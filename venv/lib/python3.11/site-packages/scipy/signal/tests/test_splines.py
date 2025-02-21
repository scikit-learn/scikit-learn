# pylint: disable=missing-docstring
import numpy as np
import pytest
from scipy._lib._array_api import xp_assert_close

from scipy.signal._spline import (
    symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd)
from scipy.signal import symiirorder1, symiirorder2


def _compute_symiirorder2_bwd_hs(k, cs, rsq, omega):
    cssq = cs * cs
    k = np.abs(k)
    rsupk = np.power(rsq, k / 2.0)

    c0 = (cssq * (1.0 + rsq) / (1.0 - rsq) /
          (1 - 2 * rsq * np.cos(2 * omega) + rsq * rsq))
    gamma = (1.0 - rsq) / (1.0 + rsq) / np.tan(omega)
    return c0 * rsupk * (np.cos(omega * k) + gamma * np.sin(omega * k))


class TestSymIIR:
    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1_ic(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Symmetrical initial conditions for a IIR filter of order 1 are:
        # x[0] + z1 * \sum{k = 0}^{n - 1} x[k] * z1^k

        # Check the initial condition for a low-pass filter
        # with coefficient b = 0.85 on a step signal. The initial condition is
        # a geometric series: 1 + b * \sum_{k = 0}^{n - 1} u[k] b^k.

        # Finding the initial condition corresponds to
        # 1. Computing the index n such that b**n < precision, which
        # corresponds to ceil(log(precision) / log(b))
        # 2. Computing the geometric series until n, this can be computed
        # using the partial sum formula: (1 - b**n) / (1 - b)
        # This holds due to the input being a step signal.
        b = 0.85
        n_exp = int(np.ceil(np.log(c_precision) / np.log(b)))
        expected = np.asarray([[(1 - b ** n_exp) / (1 - b)]], dtype=dtype)
        expected = 1 + b * expected

        # Create a step signal of size n + 1
        x = np.ones(n_exp + 1, dtype=dtype)
        xp_assert_close(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

        # Check the conditions for a exponential decreasing signal with base 2.
        # Same conditions hold, as the product of 0.5^n * 0.85^n is
        # still a geometric series
        b_d = np.asarray(b, dtype=dtype)
        expected = np.asarray(
            [[(1 - (0.5 * b_d) ** n_exp) / (1 - (0.5 * b_d))]], dtype=dtype)
        expected = 1 + b_d * expected

        # Create an exponential decreasing signal of size n + 1
        x = 2 ** -np.arange(n_exp + 1, dtype=dtype)
        xp_assert_close(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

    def test_symiir1_ic_fails(self):
        # Test that symiirorder1_ic fails whenever \sum_{n = 1}^{n} b^n > eps
        b = 0.85
        # Create a step signal of size 100
        x = np.ones(100, dtype=np.float64)

        # Compute the closed form for the geometrical series
        precision = 1 / (1 - b)
        pytest.raises(ValueError, symiirorder1_ic, x, b, precision)

        # Test that symiirorder1_ic fails when |z1| >= 1
        pytest.raises(ValueError, symiirorder1_ic, x, 1.0, -1)
        pytest.raises(ValueError, symiirorder1_ic, x, 2.0, -1)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Test for a low-pass filter with c0 = 0.15 and z1 = 0.85
        # using an unit step over 200 samples.
        c0 = 0.15
        z1 = 0.85
        n = 200
        signal = np.ones(n, dtype=dtype)

        # Find the initial condition. See test_symiir1_ic for a detailed
        # explanation
        n_exp = int(np.ceil(np.log(c_precision) / np.log(z1)))
        initial = np.asarray((1 - z1 ** n_exp) / (1 - z1), dtype=dtype)
        initial = 1 + z1 * initial

        # Forward pass
        # The transfer function for the system 1 / (1 - z1 * z^-1) when
        # applied to an unit step with initial conditions y0 is
        # 1 / (1 - z1 * z^-1) * (z^-1 / (1 - z^-1) + y0)

        # Solving the inverse Z-transform for the given expression yields:
        # y[n] = y0 * z1**n * u[n] +
        #        -z1 / (1 - z1) * z1**(k - 1) * u[k - 1] +
        #        1 / (1 - z1) * u[k - 1]
        # d is the Kronecker delta function, and u is the unit step

        # y0 * z1**n * u[n]
        pos = np.arange(n, dtype=dtype)
        comp1 = initial * z1**pos

        # -z1 / (1 - z1) * z1**(k - 1) * u[k - 1]
        comp2 = np.zeros(n, dtype=dtype)
        comp2[1:] = -z1 / (1 - z1) * z1**pos[:-1]

        # 1 / (1 - z1) * u[k - 1]
        comp3 = np.zeros(n, dtype=dtype)
        comp3[1:] = 1 / (1 - z1)

        expected_fwd = comp1 + comp2 + comp3

        # Reverse condition
        sym_cond = -c0 / (z1 - 1.0) * expected_fwd[-1]

        # Backward pass
        # The transfer function for the forward result is equivalent to
        # the forward system times c0 / (1 - z1 * z).

        # Computing a closed form for the complete expression is difficult
        # The result will be computed iteratively from the difference equation
        exp_out = np.zeros(n, dtype=dtype)
        exp_out[0] = sym_cond

        for i in range(1, n):
            exp_out[i] = c0 * expected_fwd[n - 1 - i] + z1 * exp_out[i - 1]

        exp_out = exp_out[::-1]

        out = symiirorder1(signal, c0, z1, precision)
        xp_assert_close(out, exp_out, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    def test_symiir1_values(self, dtype):
        rng = np.random.RandomState(1234)
        dtype = getattr(np, dtype)
        s = rng.uniform(size=16).astype(dtype)
        res = symiirorder1(s, 0.5, 0.1)

        # values from scipy 1.9.1
        exp_res = np.array([0.14387447, 0.35166047, 0.29735238, 0.46295986, 0.45174927,
                            0.19982875, 0.20355805, 0.47378628, 0.57232247, 0.51597393,
                           0.25935107, 0.31438554, 0.41096728, 0.4190693 , 0.25812255,
                           0.33671467], dtype=res.dtype)
        assert res.dtype == dtype
        atol = {np.float64: 1e-15, np.float32: 1e-7}[dtype]
        xp_assert_close(res, exp_res, atol=atol)

        s = s + 1j*s
        res = symiirorder1(s, 0.5, 0.1)
        assert res.dtype == np.complex64 if dtype == np.float32 else np.complex128
        xp_assert_close(res, exp_res + 1j*exp_res, atol=atol)

    @pytest.mark.parametrize(
        'dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2_initial_fwd(self, dtype, precision):
        dtype = getattr(np, dtype)
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Compute the initial conditions for a order-two symmetrical low-pass
        # filter with r = 0.5 and omega = pi / 3 for an unit step input.
        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r**2

        # The index n for the initial condition is bound from 0 to the
        # first position where sin(omega * (n + 2)) = 0 => omega * (n + 2) = pi
        # For omega = pi / 3, the maximum initial condition occurs when
        # sqrt(3) / 2 * r**n < precision.
        # => n = log(2 * sqrt(3) / 3 * precision) / log(r)
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 2
        n_exp = min(ub, lb)

        # The forward initial condition for a filter of order two is:
        # \frac{cs}{\sin(\omega)} \sum_{n = 0}^{N - 1} {
        #    r^(n + 1) \sin{\omega(n + 2)}} + cs
        # The closed expression for this sum is:
        # s[n] = 2 * r * np.cos(omega) -
        #        r**2 - r**(n + 2) * np.sin(omega * (n + 3)) / np.sin(omega) +
        #        r**(n + 3) * np.sin(omega * (n + 2)) / np.sin(omega) + cs
        fwd_initial_1 = (
            cs +
            2 * r * np.cos(omega) -
            r**2 -
            r**(n_exp + 2) * np.sin(omega * (n_exp + 3)) / np.sin(omega) +
            r**(n_exp + 3) * np.sin(omega * (n_exp + 2)) / np.sin(omega))

        # The second initial condition is given by
        # s[n] = 1 / np.sin(omega) * (
        #        r**2 * np.sin(3 * omega) -
        #        r**3 * np.sin(2 * omega) -
        #        r**(n + 3) * np.sin(omega * (n + 4)) +
        #        r**(n + 4) * np.sin(omega * (n + 3)))
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 3
        n_exp = min(ub, lb)

        fwd_initial_2 = (
            cs + cs * 2 * r * np.cos(omega) +
            (r**2 * np.sin(3 * omega) -
             r**3 * np.sin(2 * omega) -
             r**(n_exp + 3) * np.sin(omega * (n_exp + 4)) +
             r**(n_exp + 4) * np.sin(omega * (n_exp + 3))) / np.sin(omega))

        expected = np.r_[fwd_initial_1, fwd_initial_2][None, :]
        expected = expected.astype(dtype)

        n = 100
        signal = np.ones(n, dtype=dtype)

        out = symiirorder2_ic_fwd(signal, r, omega, precision)
        xp_assert_close(out, expected, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2_initial_bwd(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r * r
        a2 = 2 * r * np.cos(omega)
        a3 = -r * r

        n = 100
        signal = np.ones(n, dtype=dtype)

        # Compute initial forward conditions
        ic = symiirorder2_ic_fwd(signal, r, omega, precision)
        out = np.zeros(n + 2, dtype=dtype)
        out[:2] = ic[0]

        # Apply the forward system cs / (1 - a2 * z^-1 - a3 * z^-2))
        for i in range(2, n + 2):
            out[i] = cs * signal[i - 2] + a2 * out[i - 1] + a3 * out[i - 2]

        # Find the backward initial conditions
        ic2 = np.zeros(2, dtype=dtype)
        idx = np.arange(n)

        diff = (_compute_symiirorder2_bwd_hs(idx, cs, r * r, omega) +
                _compute_symiirorder2_bwd_hs(idx + 1, cs, r * r, omega))
        ic2_0_all = np.cumsum(diff * out[:1:-1])
        pos = np.where(diff ** 2 < c_precision)[0]
        ic2[0] = ic2_0_all[pos[0]]

        diff = (_compute_symiirorder2_bwd_hs(idx - 1, cs, r * r, omega) +
                _compute_symiirorder2_bwd_hs(idx + 2, cs, r * r, omega))
        ic2_1_all = np.cumsum(diff * out[:1:-1])
        pos = np.where(diff ** 2 < c_precision)[0]
        ic2[1] = ic2_1_all[pos[0]]

        out_ic = symiirorder2_ic_bwd(out, r, omega, precision)[0]
        xp_assert_close(out_ic, ic2, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2(self, dtype, precision):
        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r * r
        a2 = 2 * r * np.cos(omega)
        a3 = -r * r

        n = 100
        signal = np.ones(n, dtype=dtype)

        # Compute initial forward conditions
        ic = symiirorder2_ic_fwd(signal, r, omega, precision)
        out1 = np.zeros(n + 2, dtype=dtype)
        out1[:2] = ic[0]

        # Apply the forward system cs / (1 - a2 * z^-1 - a3 * z^-2))
        for i in range(2, n + 2):
            out1[i] = cs * signal[i - 2] + a2 * out1[i - 1] + a3 * out1[i - 2]

        # Find the backward initial conditions
        ic2 = symiirorder2_ic_bwd(out1, r, omega, precision)[0]

        # Apply the system cs / (1 - a2 * z - a3 * z^2)) in backwards
        exp = np.empty(n, dtype=dtype)
        exp[-2:] = ic2[::-1]

        for i in range(n - 3, -1, -1):
            exp[i] = cs * out1[i] + a2 * exp[i + 1] + a3 * exp[i + 2]

        out = symiirorder2(signal, r, omega, precision)
        xp_assert_close(out, exp, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize('dtyp', ['float32', 'float64'])
    def test_symiir2_values(self, dtyp):
        dtyp = getattr(np, dtyp)
        rng = np.random.RandomState(1234)
        s = rng.uniform(size=16).astype(dtyp)
        res = symiirorder2(s, 0.1, 0.1, precision=1e-10)

        # values from scipy 1.9.1
        exp_res = np.array([0.26572609, 0.53408018, 0.51032696, 0.72115829, 0.69486885,
           0.3649055 , 0.37349478, 0.74165032, 0.89718521, 0.80582483,
           0.46758053, 0.51898709, 0.65025605, 0.65394321, 0.45273595,
           0.53539183], dtype=dtyp)

        assert res.dtype == dtyp
        # The values in SciPy 1.14 agree with those in SciPy 1.9.1 to this
        # accuracy only. Implementation differences are twofold:
        # 1. boundary conditions are computed differently
        # 2. the filter itself uses sosfilt instead of a hardcoded iteration
        # The boundary conditions seem are tested separately (see
        # test_symiir2_initial_{fwd,bwd} above, so the difference is likely
        # due to a different way roundoff errors accumulate in the filter.
        # In that respect, sosfilt is likely doing a better job.
        xp_assert_close(res, exp_res, atol=2e-6)

        s = s + 1j*s
        with pytest.raises(TypeError):
            res = symiirorder2(s, 0.5, 0.1)

    def test_symiir1_integer_input(self):
        s = np.where(np.arange(100) % 2, -1, 1)
        expected = symiirorder1(s.astype(float), 0.5, 0.5)
        out = symiirorder1(s, 0.5, 0.5)
        xp_assert_close(out, expected)

    def test_symiir2_integer_input(self):
        s = np.where(np.arange(100) % 2, -1, 1)
        expected = symiirorder2(s.astype(float), 0.5, np.pi / 3.0)
        out = symiirorder2(s, 0.5, np.pi / 3.0)
        xp_assert_close(out, expected)
