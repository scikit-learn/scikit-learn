from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import assert_equal, assert_allclose, assert_
from scipy._lib._numpy_compat import suppress_warnings
from pytest import raises as assert_raises
import pytest

from scipy.interpolate import (BSpline, BPoly, PPoly, make_interp_spline,
        make_lsq_spline, _bspl, splev, splrep, splprep, splder, splantider,
         sproot, splint, insert)
import scipy.linalg as sl
from scipy._lib._version import NumpyVersion

from scipy.interpolate._bsplines import _not_a_knot, _augknt
import scipy.interpolate._fitpack_impl as _impl
from scipy.interpolate._fitpack import _splint


class TestBSpline(object):

    def test_ctor(self):
        # knots should be an ordered 1D array of finite real numbers
        assert_raises((TypeError, ValueError), BSpline,
                **dict(t=[1, 1.j], c=[1.], k=0))
        with np.errstate(invalid='ignore'):
            assert_raises(ValueError, BSpline, **dict(t=[1, np.nan], c=[1.], k=0))
        assert_raises(ValueError, BSpline, **dict(t=[1, np.inf], c=[1.], k=0))
        assert_raises(ValueError, BSpline, **dict(t=[1, -1], c=[1.], k=0))
        assert_raises(ValueError, BSpline, **dict(t=[[1], [1]], c=[1.], k=0))

        # for n+k+1 knots and degree k need at least n coefficients
        assert_raises(ValueError, BSpline, **dict(t=[0, 1, 2], c=[1], k=0))
        assert_raises(ValueError, BSpline,
                **dict(t=[0, 1, 2, 3, 4], c=[1., 1.], k=2))

        # non-integer orders
        assert_raises(TypeError, BSpline,
                **dict(t=[0., 0., 1., 2., 3., 4.], c=[1., 1., 1.], k="cubic"))
        assert_raises(TypeError, BSpline,
                **dict(t=[0., 0., 1., 2., 3., 4.], c=[1., 1., 1.], k=2.5))

        # basic interval cannot have measure zero (here: [1..1])
        assert_raises(ValueError, BSpline,
                **dict(t=[0., 0, 1, 1, 2, 3], c=[1., 1, 1], k=2))

        # tck vs self.tck
        n, k = 11, 3
        t = np.arange(n+k+1)
        c = np.random.random(n)
        b = BSpline(t, c, k)

        assert_allclose(t, b.t)
        assert_allclose(c, b.c)
        assert_equal(k, b.k)

    def test_tck(self):
        b = _make_random_spline()
        tck = b.tck

        assert_allclose(b.t, tck[0], atol=1e-15, rtol=1e-15)
        assert_allclose(b.c, tck[1], atol=1e-15, rtol=1e-15)
        assert_equal(b.k, tck[2])

        # b.tck is read-only
        with pytest.raises(AttributeError):
            b.tck = 'foo'

    def test_degree_0(self):
        xx = np.linspace(0, 1, 10)

        b = BSpline(t=[0, 1], c=[3.], k=0)
        assert_allclose(b(xx), 3)

        b = BSpline(t=[0, 0.35, 1], c=[3, 4], k=0)
        assert_allclose(b(xx), np.where(xx < 0.35, 3, 4))

    def test_degree_1(self):
        t = [0, 1, 2, 3, 4]
        c = [1, 2, 3]
        k = 1
        b = BSpline(t, c, k)

        x = np.linspace(1, 3, 50)
        assert_allclose(c[0]*B_012(x) + c[1]*B_012(x-1) + c[2]*B_012(x-2),
                        b(x), atol=1e-14)
        assert_allclose(splev(x, (t, c, k)), b(x), atol=1e-14)

    def test_bernstein(self):
        # a special knot vector: Bernstein polynomials
        k = 3
        t = np.asarray([0]*(k+1) + [1]*(k+1))
        c = np.asarray([1., 2., 3., 4.])
        bp = BPoly(c.reshape(-1, 1), [0, 1])
        bspl = BSpline(t, c, k)

        xx = np.linspace(-1., 2., 10)
        assert_allclose(bp(xx, extrapolate=True),
                        bspl(xx, extrapolate=True), atol=1e-14)
        assert_allclose(splev(xx, (t, c, k)),
                        bspl(xx), atol=1e-14)

    def test_rndm_naive_eval(self):
        # test random coefficient spline *on the base interval*,
        # t[k] <= x < t[-k-1]
        b = _make_random_spline()
        t, c, k = b.tck
        xx = np.linspace(t[k], t[-k-1], 50)
        y_b = b(xx)

        y_n = [_naive_eval(x, t, c, k) for x in xx]
        assert_allclose(y_b, y_n, atol=1e-14)

        y_n2 = [_naive_eval_2(x, t, c, k) for x in xx]
        assert_allclose(y_b, y_n2, atol=1e-14)

    def test_rndm_splev(self):
        b = _make_random_spline()
        t, c, k = b.tck
        xx = np.linspace(t[k], t[-k-1], 50)
        assert_allclose(b(xx), splev(xx, (t, c, k)), atol=1e-14)

    def test_rndm_splrep(self):
        np.random.seed(1234)
        x = np.sort(np.random.random(20))
        y = np.random.random(20)

        tck = splrep(x, y)
        b = BSpline(*tck)

        t, k = b.t, b.k
        xx = np.linspace(t[k], t[-k-1], 80)
        assert_allclose(b(xx), splev(xx, tck), atol=1e-14)

    def test_rndm_unity(self):
        b = _make_random_spline()
        b.c = np.ones_like(b.c)
        xx = np.linspace(b.t[b.k], b.t[-b.k-1], 100)
        assert_allclose(b(xx), 1.)

    def test_vectorization(self):
        n, k = 22, 3
        t = np.sort(np.random.random(n))
        c = np.random.random(size=(n, 6, 7))
        b = BSpline(t, c, k)
        tm, tp = t[k], t[-k-1]
        xx = tm + (tp - tm) * np.random.random((3, 4, 5))
        assert_equal(b(xx).shape, (3, 4, 5, 6, 7))

    def test_len_c(self):
        # for n+k+1 knots, only first n coefs are used.
        # and BTW this is consistent with FITPACK
        n, k = 33, 3
        t = np.sort(np.random.random(n+k+1))
        c = np.random.random(n)

        # pad coefficients with random garbage
        c_pad = np.r_[c, np.random.random(k+1)]

        b, b_pad = BSpline(t, c, k), BSpline(t, c_pad, k)

        dt = t[-1] - t[0]
        xx = np.linspace(t[0] - dt, t[-1] + dt, 50)
        assert_allclose(b(xx), b_pad(xx), atol=1e-14)
        assert_allclose(b(xx), splev(xx, (t, c, k)), atol=1e-14)
        assert_allclose(b(xx), splev(xx, (t, c_pad, k)), atol=1e-14)

    def test_endpoints(self):
        # base interval is closed
        b = _make_random_spline()
        t, _, k = b.tck
        tm, tp = t[k], t[-k-1]
        for extrap in (True, False):
            assert_allclose(b([tm, tp], extrap),
                            b([tm + 1e-10, tp - 1e-10], extrap), atol=1e-9)

    def test_continuity(self):
        # assert continuity at internal knots
        b = _make_random_spline()
        t, _, k = b.tck
        assert_allclose(b(t[k+1:-k-1] - 1e-10), b(t[k+1:-k-1] + 1e-10),
                atol=1e-9)

    def test_extrap(self):
        b = _make_random_spline()
        t, c, k = b.tck
        dt = t[-1] - t[0]
        xx = np.linspace(t[k] - dt, t[-k-1] + dt, 50)
        mask = (t[k] < xx) & (xx < t[-k-1])

        # extrap has no effect within the base interval
        assert_allclose(b(xx[mask], extrapolate=True),
                        b(xx[mask], extrapolate=False))

        # extrapolated values agree with FITPACK
        assert_allclose(b(xx, extrapolate=True),
                splev(xx, (t, c, k), ext=0))

    def test_default_extrap(self):
        # BSpline defaults to extrapolate=True
        b = _make_random_spline()
        t, _, k = b.tck
        xx = [t[0] - 1, t[-1] + 1]
        yy = b(xx)
        assert_(not np.all(np.isnan(yy)))

    def test_periodic_extrap(self):
        np.random.seed(1234)
        t = np.sort(np.random.random(8))
        c = np.random.random(4)
        k = 3
        b = BSpline(t, c, k, extrapolate='periodic')
        n = t.size - (k + 1)

        dt = t[-1] - t[0]
        xx = np.linspace(t[k] - dt, t[n] + dt, 50)
        xy = t[k] + (xx - t[k]) % (t[n] - t[k])
        assert_allclose(b(xx), splev(xy, (t, c, k)))

        # Direct check
        xx = [-1, 0, 0.5, 1]
        xy = t[k] + (xx - t[k]) % (t[n] - t[k])
        assert_equal(b(xx, extrapolate='periodic'), b(xy, extrapolate=True))

    def test_ppoly(self):
        b = _make_random_spline()
        t, c, k = b.tck
        pp = PPoly.from_spline((t, c, k))

        xx = np.linspace(t[k], t[-k], 100)
        assert_allclose(b(xx), pp(xx), atol=1e-14, rtol=1e-14)

    def test_derivative_rndm(self):
        b = _make_random_spline()
        t, c, k = b.tck
        xx = np.linspace(t[0], t[-1], 50)
        xx = np.r_[xx, t]

        for der in range(1, k+1):
            yd = splev(xx, (t, c, k), der=der)
            assert_allclose(yd, b(xx, nu=der), atol=1e-14)

        # higher derivatives all vanish
        assert_allclose(b(xx, nu=k+1), 0, atol=1e-14)

    def test_derivative_jumps(self):
        # example from de Boor, Chap IX, example (24)
        # NB: knots augmented & corresp coefs are zeroed out
        # in agreement with the convention (29)
        k = 2
        t = [-1, -1, 0, 1, 1, 3, 4, 6, 6, 6, 7, 7]
        np.random.seed(1234)
        c = np.r_[0, 0, np.random.random(5), 0, 0]
        b = BSpline(t, c, k)

        # b is continuous at x != 6 (triple knot)
        x = np.asarray([1, 3, 4, 6])
        assert_allclose(b(x[x != 6] - 1e-10),
                        b(x[x != 6] + 1e-10))
        assert_(not np.allclose(b(6.-1e-10), b(6+1e-10)))

        # 1st derivative jumps at double knots, 1 & 6:
        x0 = np.asarray([3, 4])
        assert_allclose(b(x0 - 1e-10, nu=1),
                        b(x0 + 1e-10, nu=1))
        x1 = np.asarray([1, 6])
        assert_(not np.all(np.allclose(b(x1 - 1e-10, nu=1),
                                       b(x1 + 1e-10, nu=1))))

        # 2nd derivative is not guaranteed to be continuous either
        assert_(not np.all(np.allclose(b(x - 1e-10, nu=2),
                                       b(x + 1e-10, nu=2))))

    def test_basis_element_quadratic(self):
        xx = np.linspace(-1, 4, 20)
        b = BSpline.basis_element(t=[0, 1, 2, 3])
        assert_allclose(b(xx),
                        splev(xx, (b.t, b.c, b.k)), atol=1e-14)
        assert_allclose(b(xx),
                        B_0123(xx), atol=1e-14)

        b = BSpline.basis_element(t=[0, 1, 1, 2])
        xx = np.linspace(0, 2, 10)
        assert_allclose(b(xx),
                np.where(xx < 1, xx*xx, (2.-xx)**2), atol=1e-14)

    def test_basis_element_rndm(self):
        b = _make_random_spline()
        t, c, k = b.tck
        xx = np.linspace(t[k], t[-k-1], 20)
        assert_allclose(b(xx), _sum_basis_elements(xx, t, c, k), atol=1e-14)

    def test_cmplx(self):
        b = _make_random_spline()
        t, c, k = b.tck
        cc = c * (1. + 3.j)

        b = BSpline(t, cc, k)
        b_re = BSpline(t, b.c.real, k)
        b_im = BSpline(t, b.c.imag, k)

        xx = np.linspace(t[k], t[-k-1], 20)
        assert_allclose(b(xx).real, b_re(xx), atol=1e-14)
        assert_allclose(b(xx).imag, b_im(xx), atol=1e-14)

    def test_nan(self):
        # nan in, nan out.
        b = BSpline.basis_element([0, 1, 1, 2])
        assert_(np.isnan(b(np.nan)))

    def test_derivative_method(self):
        b = _make_random_spline(k=5)
        t, c, k = b.tck
        b0 = BSpline(t, c, k)
        xx = np.linspace(t[k], t[-k-1], 20)
        for j in range(1, k):
            b = b.derivative()
            assert_allclose(b0(xx, j), b(xx), atol=1e-12, rtol=1e-12)

    def test_antiderivative_method(self):
        b = _make_random_spline()
        t, c, k = b.tck
        xx = np.linspace(t[k], t[-k-1], 20)
        assert_allclose(b.antiderivative().derivative()(xx),
                        b(xx), atol=1e-14, rtol=1e-14)

        # repeat with n-D array for c
        c = np.c_[c, c, c]
        c = np.dstack((c, c))
        b = BSpline(t, c, k)
        assert_allclose(b.antiderivative().derivative()(xx),
                        b(xx), atol=1e-14, rtol=1e-14)

    def test_integral(self):
        b = BSpline.basis_element([0, 1, 2])  # x for x < 1 else 2 - x
        assert_allclose(b.integrate(0, 1), 0.5)
        assert_allclose(b.integrate(1, 0), -1 * 0.5)
        assert_allclose(b.integrate(1, 0), -0.5)

        # extrapolate or zeros outside of [0, 2]; default is yes
        assert_allclose(b.integrate(-1, 1), 0)
        assert_allclose(b.integrate(-1, 1, extrapolate=True), 0)
        assert_allclose(b.integrate(-1, 1, extrapolate=False), 0.5)
        assert_allclose(b.integrate(1, -1, extrapolate=False), -1 * 0.5)

        # Test ``_fitpack._splint()``
        t, c, k = b.tck
        assert_allclose(b.integrate(1, -1, extrapolate=False),
                        _splint(t, c, k, 1, -1)[0])

        # Test ``extrapolate='periodic'``.
        b.extrapolate = 'periodic'
        i = b.antiderivative()
        period_int = i(2) - i(0)

        assert_allclose(b.integrate(0, 2), period_int)
        assert_allclose(b.integrate(2, 0), -1 * period_int)
        assert_allclose(b.integrate(-9, -7), period_int)
        assert_allclose(b.integrate(-8, -4), 2 * period_int)

        assert_allclose(b.integrate(0.5, 1.5), i(1.5) - i(0.5))
        assert_allclose(b.integrate(1.5, 3), i(1) - i(0) + i(2) - i(1.5))
        assert_allclose(b.integrate(1.5 + 12, 3 + 12),
                        i(1) - i(0) + i(2) - i(1.5))
        assert_allclose(b.integrate(1.5, 3 + 12),
                        i(1) - i(0) + i(2) - i(1.5) + 6 * period_int)

        assert_allclose(b.integrate(0, -1), i(0) - i(1))
        assert_allclose(b.integrate(-9, -10), i(0) - i(1))
        assert_allclose(b.integrate(0, -9), i(1) - i(2) - 4 * period_int)

    def test_integrate_ppoly(self):
        # test .integrate method to be consistent with PPoly.integrate
        x = [0, 1, 2, 3, 4]
        b = make_interp_spline(x, x)
        b.extrapolate = 'periodic'
        p = PPoly.from_spline(b)

        for x0, x1 in [(-5, 0.5), (0.5, 5), (-4, 13)]:
            assert_allclose(b.integrate(x0, x1),
                            p.integrate(x0, x1))

    def test_subclassing(self):
        # classmethods should not decay to the base class
        class B(BSpline):
            pass

        b = B.basis_element([0, 1, 2, 2])
        assert_equal(b.__class__, B)
        assert_equal(b.derivative().__class__, B)
        assert_equal(b.antiderivative().__class__, B)

    def test_axis(self):
        n, k = 22, 3
        t = np.linspace(0, 1, n + k + 1)
        sh0 = [6, 7, 8]
        for axis in range(4):
            sh = sh0[:]
            sh.insert(axis, n)   # [22, 6, 7, 8] etc
            c = np.random.random(size=sh)
            b = BSpline(t, c, k, axis=axis)
            assert_equal(b.c.shape,
                         [sh[axis],] + sh[:axis] + sh[axis+1:])

            xp = np.random.random((3, 4, 5))
            assert_equal(b(xp).shape,
                         sh[:axis] + list(xp.shape) + sh[axis+1:])

            #0 <= axis < c.ndim
            for ax in [-1, c.ndim]:
                assert_raises(ValueError, BSpline, **dict(t=t, c=c, k=k, axis=ax))

            # derivative, antiderivative keeps the axis
            for b1 in [BSpline(t, c, k, axis=axis).derivative(),
                       BSpline(t, c, k, axis=axis).derivative(2),
                       BSpline(t, c, k, axis=axis).antiderivative(),
                       BSpline(t, c, k, axis=axis).antiderivative(2)]:
                assert_equal(b1.axis, b.axis)


def test_knots_multiplicity():
    # Take a spline w/ random coefficients, throw in knots of varying
    # multiplicity.

    def check_splev(b, j, der=0, atol=1e-14, rtol=1e-14):
        # check evaluations against FITPACK, incl extrapolations
        t, c, k = b.tck
        x = np.unique(t)
        x = np.r_[t[0]-0.1, 0.5*(x[1:] + x[:1]), t[-1]+0.1]
        assert_allclose(splev(x, (t, c, k), der), b(x, der),
                atol=atol, rtol=rtol, err_msg='der = %s  k = %s' % (der, b.k))

    # test loop itself
    # [the index `j` is for interpreting the traceback in case of a failure]
    for k in [1, 2, 3, 4, 5]:
        b = _make_random_spline(k=k)
        for j, b1 in enumerate(_make_multiples(b)):
            check_splev(b1, j)
            for der in range(1, k+1):
                check_splev(b1, j, der, 1e-12, 1e-12)


### stolen from @pv, verbatim
def _naive_B(x, k, i, t):
    """
    Naive way to compute B-spline basis functions. Useful only for testing!
    computes B(x; t[i],..., t[i+k+1])
    """
    if k == 0:
        return 1.0 if t[i] <= x < t[i+1] else 0.0
    if t[i+k] == t[i]:
        c1 = 0.0
    else:
        c1 = (x - t[i])/(t[i+k] - t[i]) * _naive_B(x, k-1, i, t)
    if t[i+k+1] == t[i+1]:
        c2 = 0.0
    else:
        c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * _naive_B(x, k-1, i+1, t)
    return (c1 + c2)


### stolen from @pv, verbatim
def _naive_eval(x, t, c, k):
    """
    Naive B-spline evaluation. Useful only for testing!
    """
    if x == t[k]:
        i = k
    else:
        i = np.searchsorted(t, x) - 1
    assert t[i] <= x <= t[i+1]
    assert i >= k and i < len(t) - k
    return sum(c[i-j] * _naive_B(x, k, i-j, t) for j in range(0, k+1))


def _naive_eval_2(x, t, c, k):
    """Naive B-spline evaluation, another way."""
    n = len(t) - (k+1)
    assert n >= k+1
    assert len(c) >= n
    assert t[k] <= x <= t[n]
    return sum(c[i] * _naive_B(x, k, i, t) for i in range(n))


def _sum_basis_elements(x, t, c, k):
    n = len(t) - (k+1)
    assert n >= k+1
    assert len(c) >= n
    s = 0.
    for i in range(n):
        b = BSpline.basis_element(t[i:i+k+2], extrapolate=False)(x)
        s += c[i] * np.nan_to_num(b)   # zero out out-of-bounds elements
    return s


def B_012(x):
    """ A linear B-spline function B(x | 0, 1, 2)."""
    x = np.atleast_1d(x)
    return np.piecewise(x, [(x < 0) | (x > 2),
                            (x >= 0) & (x < 1),
                            (x >= 1) & (x <= 2)],
                           [lambda x: 0., lambda x: x, lambda x: 2.-x])


def B_0123(x, der=0):
    """A quadratic B-spline function B(x | 0, 1, 2, 3)."""
    x = np.atleast_1d(x)
    conds = [x < 1, (x > 1) & (x < 2), x > 2]
    if der == 0:
        funcs = [lambda x: x*x/2.,
                 lambda x: 3./4 - (x-3./2)**2,
                 lambda x: (3.-x)**2 / 2]
    elif der == 2:
        funcs = [lambda x: 1.,
                 lambda x: -2.,
                 lambda x: 1.]
    else:
        raise ValueError('never be here: der=%s' % der)
    pieces = np.piecewise(x, conds, funcs)
    return pieces


def _make_random_spline(n=35, k=3):
    np.random.seed(123)
    t = np.sort(np.random.random(n+k+1))
    c = np.random.random(n)
    return BSpline.construct_fast(t, c, k)


def _make_multiples(b):
    """Increase knot multiplicity."""
    c, k = b.c, b.k

    t1 = b.t.copy()
    t1[17:19] = t1[17]
    t1[22] = t1[21]
    yield BSpline(t1, c, k)

    t1 = b.t.copy()
    t1[:k+1] = t1[0]
    yield BSpline(t1, c, k)

    t1 = b.t.copy()
    t1[-k-1:] = t1[-1]
    yield BSpline(t1, c, k)


class TestInterop(object):
    #
    # Test that FITPACK-based spl* functions can deal with BSpline objects
    #
    def setup_method(self):
        xx = np.linspace(0, 4.*np.pi, 41)
        yy = np.cos(xx)
        b = make_interp_spline(xx, yy)
        self.tck = (b.t, b.c, b.k)
        self.xx, self.yy, self.b = xx, yy, b

        self.xnew = np.linspace(0, 4.*np.pi, 21)

        c2 = np.c_[b.c, b.c, b.c]
        self.c2 = np.dstack((c2, c2))
        self.b2 = BSpline(b.t, self.c2, b.k)

    def test_splev(self):
        xnew, b, b2 = self.xnew, self.b, self.b2

        # check that splev works with 1D array of coefficients
        # for array and scalar `x`
        assert_allclose(splev(xnew, b),
                        b(xnew), atol=1e-15, rtol=1e-15)
        assert_allclose(splev(xnew, b.tck),
                        b(xnew), atol=1e-15, rtol=1e-15)
        assert_allclose([splev(x, b) for x in xnew],
                        b(xnew), atol=1e-15, rtol=1e-15)

        # With n-D coefficients, there's a quirck:
        # splev(x, BSpline) is equivalent to BSpline(x)
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning,
                       "Calling splev.. with BSpline objects with c.ndim > 1 is not recommended.")
            assert_allclose(splev(xnew, b2), b2(xnew), atol=1e-15, rtol=1e-15)

        # However, splev(x, BSpline.tck) needs some transposes. This is because
        # BSpline interpolates along the first axis, while the legacy FITPACK
        # wrapper does list(map(...)) which effectively interpolates along the
        # last axis. Like so:
        sh = tuple(range(1, b2.c.ndim)) + (0,)   # sh = (1, 2, 0)
        cc = b2.c.transpose(sh)
        tck = (b2.t, cc, b2.k)
        assert_allclose(splev(xnew, tck),
                        b2(xnew).transpose(sh), atol=1e-15, rtol=1e-15)

    def test_splrep(self):
        x, y = self.xx, self.yy
        # test that "new" splrep is equivalent to _impl.splrep
        tck = splrep(x, y)
        t, c, k = _impl.splrep(x, y)
        assert_allclose(tck[0], t, atol=1e-15)
        assert_allclose(tck[1], c, atol=1e-15)
        assert_equal(tck[2], k)

        # also cover the `full_output=True` branch
        tck_f, _, _, _ = splrep(x, y, full_output=True)
        assert_allclose(tck_f[0], t, atol=1e-15)
        assert_allclose(tck_f[1], c, atol=1e-15)
        assert_equal(tck_f[2], k)

        # test that the result of splrep roundtrips with splev:
        # evaluate the spline on the original `x` points
        yy = splev(x, tck)
        assert_allclose(y, yy, atol=1e-15)

        # ... and also it roundtrips if wrapped in a BSpline
        b = BSpline(*tck)
        assert_allclose(y, b(x), atol=1e-15)

    @pytest.mark.xfail(NumpyVersion(np.__version__) < '1.14.0',
                       reason='requires NumPy >= 1.14.0')
    def test_splrep_errors(self):
        # test that both "old" and "new" splrep raise for an n-D ``y`` array
        # with n > 1
        x, y = self.xx, self.yy
        y2 = np.c_[y, y]
        with assert_raises(ValueError):
            splrep(x, y2)
        with assert_raises(ValueError):
            _impl.splrep(x, y2)

        # input below minimum size
        with assert_raises(TypeError, match="m > k must hold"):
            splrep(x[:3], y[:3])
        with assert_raises(TypeError, match="m > k must hold"):
            _impl.splrep(x[:3], y[:3])

    def test_splprep(self):
        x = np.arange(15).reshape((3, 5))
        b, u = splprep(x)
        tck, u1 = _impl.splprep(x)

        # test the roundtrip with splev for both "old" and "new" output
        assert_allclose(u, u1, atol=1e-15)
        assert_allclose(splev(u, b), x, atol=1e-15)
        assert_allclose(splev(u, tck), x, atol=1e-15)

        # cover the ``full_output=True`` branch
        (b_f, u_f), _, _, _ = splprep(x, s=0, full_output=True)
        assert_allclose(u, u_f, atol=1e-15)
        assert_allclose(splev(u_f, b_f), x, atol=1e-15)

    def test_splprep_errors(self):
        # test that both "old" and "new" code paths raise for x.ndim > 2
        x = np.arange(3*4*5).reshape((3, 4, 5))
        with assert_raises(ValueError, match="too many values to unpack"):
            splprep(x)
        with assert_raises(ValueError, match="too many values to unpack"):
            _impl.splprep(x)

        # input below minimum size
        x = np.linspace(0, 40, num=3)
        with assert_raises(TypeError, match="m > k must hold"):
            splprep([x])
        with assert_raises(TypeError, match="m > k must hold"):
            _impl.splprep([x])

        # automatically calculated parameters are non-increasing
        # see gh-7589
        x = [-50.49072266, -50.49072266, -54.49072266, -54.49072266]
        with assert_raises(ValueError, match="Invalid inputs"):
            splprep([x])
        with assert_raises(ValueError, match="Invalid inputs"):
            _impl.splprep([x])

        # given non-increasing parameter values u
        x = [1, 3, 2, 4]
        u = [0, 0.3, 0.2, 1]
        with assert_raises(ValueError, match="Invalid inputs"):
            splprep(*[[x], None, u])

    def test_sproot(self):
        b, b2 = self.b, self.b2
        roots = np.array([0.5, 1.5, 2.5, 3.5])*np.pi
        # sproot accepts a BSpline obj w/ 1D coef array
        assert_allclose(sproot(b), roots, atol=1e-7, rtol=1e-7)
        assert_allclose(sproot((b.t, b.c, b.k)), roots, atol=1e-7, rtol=1e-7)

        # ... and deals with trailing dimensions if coef array is n-D
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning,
                       "Calling sproot.. with BSpline objects with c.ndim > 1 is not recommended.")
            r = sproot(b2, mest=50)
        r = np.asarray(r)

        assert_equal(r.shape, (3, 2, 4))
        assert_allclose(r - roots, 0, atol=1e-12)

        # and legacy behavior is preserved for a tck tuple w/ n-D coef
        c2r = b2.c.transpose(1, 2, 0)
        rr = np.asarray(sproot((b2.t, c2r, b2.k), mest=50))
        assert_equal(rr.shape, (3, 2, 4))
        assert_allclose(rr - roots, 0, atol=1e-12)

    def test_splint(self):
        # test that splint accepts BSpline objects
        b, b2 = self.b, self.b2
        assert_allclose(splint(0, 1, b),
                        splint(0, 1, b.tck), atol=1e-14)
        assert_allclose(splint(0, 1, b),
                        b.integrate(0, 1), atol=1e-14)

        # ... and deals with n-D arrays of coefficients
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning,
                       "Calling splint.. with BSpline objects with c.ndim > 1 is not recommended.")
            assert_allclose(splint(0, 1, b2), b2.integrate(0, 1), atol=1e-14)

        # and the legacy behavior is preserved for a tck tuple w/ n-D coef
        c2r = b2.c.transpose(1, 2, 0)
        integr = np.asarray(splint(0, 1, (b2.t, c2r, b2.k)))
        assert_equal(integr.shape, (3, 2))
        assert_allclose(integr,
                        splint(0, 1, b), atol=1e-14)

    def test_splder(self):
        for b in [self.b, self.b2]:
            # pad the c array (FITPACK convention)
            ct = len(b.t) - len(b.c)
            if ct > 0:
                b.c = np.r_[b.c, np.zeros((ct,) + b.c.shape[1:])]

            for n in [1, 2, 3]:
                bd = splder(b)
                tck_d = _impl.splder((b.t, b.c, b.k))
                assert_allclose(bd.t, tck_d[0], atol=1e-15)
                assert_allclose(bd.c, tck_d[1], atol=1e-15)
                assert_equal(bd.k, tck_d[2])
                assert_(isinstance(bd, BSpline))
                assert_(isinstance(tck_d, tuple))  # back-compat: tck in and out

    def test_splantider(self):
        for b in [self.b, self.b2]:
            # pad the c array (FITPACK convention)
            ct = len(b.t) - len(b.c)
            if ct > 0:
                b.c = np.r_[b.c, np.zeros((ct,) + b.c.shape[1:])]

            for n in [1, 2, 3]:
                bd = splantider(b)
                tck_d = _impl.splantider((b.t, b.c, b.k))
                assert_allclose(bd.t, tck_d[0], atol=1e-15)
                assert_allclose(bd.c, tck_d[1], atol=1e-15)
                assert_equal(bd.k, tck_d[2])
                assert_(isinstance(bd, BSpline))
                assert_(isinstance(tck_d, tuple))  # back-compat: tck in and out

    def test_insert(self):
        b, b2, xx = self.b, self.b2, self.xx

        j = b.t.size // 2
        tn = 0.5*(b.t[j] + b.t[j+1])

        bn, tck_n = insert(tn, b), insert(tn, (b.t, b.c, b.k))
        assert_allclose(splev(xx, bn),
                        splev(xx, tck_n), atol=1e-15)
        assert_(isinstance(bn, BSpline))
        assert_(isinstance(tck_n, tuple))   # back-compat: tck in, tck out

        # for n-D array of coefficients, BSpline.c needs to be transposed
        # after that, the results are equivalent.
        sh = tuple(range(b2.c.ndim))
        c_ = b2.c.transpose(sh[1:] + (0,))
        tck_n2 = insert(tn, (b2.t, c_, b2.k))

        bn2 = insert(tn, b2)

        # need a transpose for comparing the results, cf test_splev
        assert_allclose(np.asarray(splev(xx, tck_n2)).transpose(2, 0, 1),
                        bn2(xx), atol=1e-15)
        assert_(isinstance(bn2, BSpline))
        assert_(isinstance(tck_n2, tuple))   # back-compat: tck in, tck out


class TestInterp(object):
    #
    # Test basic ways of constructing interpolating splines.
    #
    xx = np.linspace(0., 2.*np.pi)
    yy = np.sin(xx)

    def test_non_int_order(self):
        with assert_raises(TypeError):
            make_interp_spline(self.xx, self.yy, k=2.5)

    def test_order_0(self):
        b = make_interp_spline(self.xx, self.yy, k=0)
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        b = make_interp_spline(self.xx, self.yy, k=0, axis=-1)
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)

    def test_linear(self):
        b = make_interp_spline(self.xx, self.yy, k=1)
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        b = make_interp_spline(self.xx, self.yy, k=1, axis=-1)
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)

    def test_not_a_knot(self):
        for k in [3, 5]:
            b = make_interp_spline(self.xx, self.yy, k)
            assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)

    def test_quadratic_deriv(self):
        der = [(1, 8.)]  # order, value: f'(x) = 8.

        # derivative at right-hand edge
        b = make_interp_spline(self.xx, self.yy, k=2, bc_type=(None, der))
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        assert_allclose(b(self.xx[-1], 1), der[0][1], atol=1e-14, rtol=1e-14)

        # derivative at left-hand edge
        b = make_interp_spline(self.xx, self.yy, k=2, bc_type=(der, None))
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        assert_allclose(b(self.xx[0], 1), der[0][1], atol=1e-14, rtol=1e-14)

    def test_cubic_deriv(self):
        k = 3

        # first derivatives at left & right edges:
        der_l, der_r = [(1, 3.)], [(1, 4.)]
        b = make_interp_spline(self.xx, self.yy, k, bc_type=(der_l, der_r))
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        assert_allclose([b(self.xx[0], 1), b(self.xx[-1], 1)],
                        [der_l[0][1], der_r[0][1]], atol=1e-14, rtol=1e-14)

        # 'natural' cubic spline, zero out 2nd derivatives at the boundaries
        der_l, der_r = [(2, 0)], [(2, 0)]
        b = make_interp_spline(self.xx, self.yy, k, bc_type=(der_l, der_r))
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)

    def test_quintic_derivs(self):
        k, n = 5, 7
        x = np.arange(n).astype(np.float_)
        y = np.sin(x)
        der_l = [(1, -12.), (2, 1)]
        der_r = [(1, 8.), (2, 3.)]
        b = make_interp_spline(x, y, k=k, bc_type=(der_l, der_r))
        assert_allclose(b(x), y, atol=1e-14, rtol=1e-14)
        assert_allclose([b(x[0], 1), b(x[0], 2)],
                        [val for (nu, val) in der_l])
        assert_allclose([b(x[-1], 1), b(x[-1], 2)],
                        [val for (nu, val) in der_r])

    @pytest.mark.xfail(reason='unstable')
    def test_cubic_deriv_unstable(self):
        # 1st and 2nd derivative at x[0], no derivative information at x[-1]
        # The problem is not that it fails [who would use this anyway],
        # the problem is that it fails *silently*, and I've no idea
        # how to detect this sort of instability.
        # In this particular case: it's OK for len(t) < 20, goes haywire
        # at larger `len(t)`.
        k = 3
        t = _augknt(self.xx, k)

        der_l = [(1, 3.), (2, 4.)]
        b = make_interp_spline(self.xx, self.yy, k, t, bc_type=(der_l, None))
        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)

    def test_knots_not_data_sites(self):
        # Knots need not coincide with the data sites.
        # use a quadratic spline, knots are at data averages,
        # two additional constraints are zero 2nd derivs at edges
        k = 2
        t = np.r_[(self.xx[0],)*(k+1),
                  (self.xx[1:] + self.xx[:-1]) / 2.,
                  (self.xx[-1],)*(k+1)]
        b = make_interp_spline(self.xx, self.yy, k, t,
                               bc_type=([(2, 0)], [(2, 0)]))

        assert_allclose(b(self.xx), self.yy, atol=1e-14, rtol=1e-14)
        assert_allclose([b(self.xx[0], 2), b(self.xx[-1], 2)], [0., 0.],
                atol=1e-14)

    def test_minimum_points_and_deriv(self):
        # interpolation of f(x) = x**3 between 0 and 1. f'(x) = 3 * xx**2 and
        # f'(0) = 0, f'(1) = 3.
        k = 3
        x = [0., 1.]
        y = [0., 1.]
        b = make_interp_spline(x, y, k, bc_type=([(1, 0.)], [(1, 3.)]))

        xx = np.linspace(0., 1.)
        yy = xx**3
        assert_allclose(b(xx), yy, atol=1e-14, rtol=1e-14)

    def test_deriv_spec(self):
        # If one of the derivatives is omitted, the spline definition is
        # incomplete.
        x = y = [1.0, 2, 3, 4, 5, 6]

        with assert_raises(ValueError):
            make_interp_spline(x, y, bc_type=([(1, 0.)], None))

        with assert_raises(ValueError):
            make_interp_spline(x, y, bc_type=(1, 0.))

        with assert_raises(ValueError):
            make_interp_spline(x, y, bc_type=[(1, 0.)])

        with assert_raises(ValueError):
            make_interp_spline(x, y, bc_type=42)

        # CubicSpline expects`bc_type=(left_pair, right_pair)`, while
        # here we expect `bc_type=(iterable, iterable)`.
        l, r = (1, 0.0), (1, 0.0)
        with assert_raises(ValueError):
            make_interp_spline(x, y, bc_type=(l, r))

    def test_complex(self):
        k = 3
        xx = self.xx
        yy = self.yy + 1.j*self.yy

        # first derivatives at left & right edges:
        der_l, der_r = [(1, 3.j)], [(1, 4.+2.j)]
        b = make_interp_spline(xx, yy, k, bc_type=(der_l, der_r))
        assert_allclose(b(xx), yy, atol=1e-14, rtol=1e-14)
        assert_allclose([b(xx[0], 1), b(xx[-1], 1)],
                        [der_l[0][1], der_r[0][1]], atol=1e-14, rtol=1e-14)

        # also test zero and first order
        for k in (0, 1):
            b = make_interp_spline(xx, yy, k=k)
            assert_allclose(b(xx), yy, atol=1e-14, rtol=1e-14)

    def test_int_xy(self):
        x = np.arange(10).astype(np.int_)
        y = np.arange(10).astype(np.int_)

        # cython chokes on "buffer type mismatch" (construction) or
        # "no matching signature found" (evaluation)
        for k in (0, 1, 2, 3):
            b = make_interp_spline(x, y, k=k)
            b(x)

    def test_sliced_input(self):
        # cython code chokes on non C contiguous arrays
        xx = np.linspace(-1, 1, 100)

        x = xx[::5]
        y = xx[::5]

        for k in (0, 1, 2, 3):
            make_interp_spline(x, y, k=k)

    def test_check_finite(self):
        # check_finite defaults to True; nans and such trigger a ValueError
        x = np.arange(10).astype(float)
        y = x**2

        for z in [np.nan, np.inf, -np.inf]:
            y[-1] = z
            assert_raises(ValueError, make_interp_spline, x, y)

    @pytest.mark.parametrize('k', [1, 2, 3, 5])
    def test_list_input(self, k):
        # regression test for gh-8714: TypeError for x, y being lists and k=2
        x = list(range(10))
        y = [a**2 for a in x]
        make_interp_spline(x, y, k=k)

    def test_multiple_rhs(self):
        yy = np.c_[np.sin(self.xx), np.cos(self.xx)]
        der_l = [(1, [1., 2.])]
        der_r = [(1, [3., 4.])]

        b = make_interp_spline(self.xx, yy, k=3, bc_type=(der_l, der_r))
        assert_allclose(b(self.xx), yy, atol=1e-14, rtol=1e-14)
        assert_allclose(b(self.xx[0], 1), der_l[0][1], atol=1e-14, rtol=1e-14)
        assert_allclose(b(self.xx[-1], 1), der_r[0][1], atol=1e-14, rtol=1e-14)

    def test_shapes(self):
        np.random.seed(1234)
        k, n = 3, 22
        x = np.sort(np.random.random(size=n))
        y = np.random.random(size=(n, 5, 6, 7))

        b = make_interp_spline(x, y, k)
        assert_equal(b.c.shape, (n, 5, 6, 7))

        # now throw in some derivatives
        d_l = [(1, np.random.random((5, 6, 7)))]
        d_r = [(1, np.random.random((5, 6, 7)))]
        b = make_interp_spline(x, y, k, bc_type=(d_l, d_r))
        assert_equal(b.c.shape, (n + k - 1, 5, 6, 7))

    def test_string_aliases(self):
        yy = np.sin(self.xx)

        # a single string is duplicated
        b1 = make_interp_spline(self.xx, yy, k=3, bc_type='natural')
        b2 = make_interp_spline(self.xx, yy, k=3, bc_type=([(2, 0)], [(2, 0)]))
        assert_allclose(b1.c, b2.c, atol=1e-15)

        # two strings are handled
        b1 = make_interp_spline(self.xx, yy, k=3,
                                bc_type=('natural', 'clamped'))
        b2 = make_interp_spline(self.xx, yy, k=3,
                                bc_type=([(2, 0)], [(1, 0)]))
        assert_allclose(b1.c, b2.c, atol=1e-15)

        # one-sided BCs are OK
        b1 = make_interp_spline(self.xx, yy, k=2, bc_type=(None, 'clamped'))
        b2 = make_interp_spline(self.xx, yy, k=2, bc_type=(None, [(1, 0.0)]))
        assert_allclose(b1.c, b2.c, atol=1e-15)

        # 'not-a-knot' is equivalent to None
        b1 = make_interp_spline(self.xx, yy, k=3, bc_type='not-a-knot')
        b2 = make_interp_spline(self.xx, yy, k=3, bc_type=None)
        assert_allclose(b1.c, b2.c, atol=1e-15)

        # unknown strings do not pass
        with assert_raises(ValueError):
            make_interp_spline(self.xx, yy, k=3, bc_type='typo')

        # string aliases are handled for 2D values
        yy = np.c_[np.sin(self.xx), np.cos(self.xx)]
        der_l = [(1, [0., 0.])]
        der_r = [(2, [0., 0.])]
        b2 = make_interp_spline(self.xx, yy, k=3, bc_type=(der_l, der_r))
        b1 = make_interp_spline(self.xx, yy, k=3,
                                bc_type=('clamped', 'natural'))
        assert_allclose(b1.c, b2.c, atol=1e-15)

        # ... and for n-D values:
        np.random.seed(1234)
        k, n = 3, 22
        x = np.sort(np.random.random(size=n))
        y = np.random.random(size=(n, 5, 6, 7))

        # now throw in some derivatives
        d_l = [(1, np.zeros((5, 6, 7)))]
        d_r = [(1, np.zeros((5, 6, 7)))]
        b1 = make_interp_spline(x, y, k, bc_type=(d_l, d_r))
        b2 = make_interp_spline(x, y, k, bc_type='clamped')
        assert_allclose(b1.c, b2.c, atol=1e-15)

    def test_full_matrix(self):
        np.random.seed(1234)
        k, n = 3, 7
        x = np.sort(np.random.random(size=n))
        y = np.random.random(size=n)
        t = _not_a_knot(x, k)

        b = make_interp_spline(x, y, k, t)
        cf = make_interp_full_matr(x, y, t, k)
        assert_allclose(b.c, cf, atol=1e-14, rtol=1e-14)


def make_interp_full_matr(x, y, t, k):
    """Assemble an spline order k with knots t to interpolate
    y(x) using full matrices.
    Not-a-knot BC only.

    This routine is here for testing only (even though it's functional).
    """
    assert x.size == y.size
    assert t.size == x.size + k + 1
    n = x.size

    A = np.zeros((n, n), dtype=np.float_)

    for j in range(n):
        xval = x[j]
        if xval == t[k]:
            left = k
        else:
            left = np.searchsorted(t, xval) - 1

        # fill a row
        bb = _bspl.evaluate_all_bspl(t, k, xval, left)
        A[j, left-k:left+1] = bb

    c = sl.solve(A, y)
    return c


### XXX: 'periodic' interp spline using full matrices
def make_interp_per_full_matr(x, y, t, k):
    x, y, t = map(np.asarray, (x, y, t))

    n = x.size
    nt = t.size - k - 1

    # have `n` conditions for `nt` coefficients; need nt-n derivatives
    assert nt - n == k - 1

    # LHS: the collocation matrix + derivatives at edges
    A = np.zeros((nt, nt), dtype=np.float_)

    # derivatives at x[0]:
    offset = 0

    if x[0] == t[k]:
        left = k
    else:
        left = np.searchsorted(t, x[0]) - 1

    if x[-1] == t[k]:
        left2 = k
    else:
        left2 = np.searchsorted(t, x[-1]) - 1

    for i in range(k-1):
        bb = _bspl.evaluate_all_bspl(t, k, x[0], left, nu=i+1)
        A[i, left-k:left+1] = bb
        bb = _bspl.evaluate_all_bspl(t, k, x[-1], left2, nu=i+1)
        A[i, left2-k:left2+1] = -bb
        offset += 1

    # RHS
    y = np.r_[[0]*(k-1), y]

    # collocation matrix
    for j in range(n):
        xval = x[j]
        # find interval
        if xval == t[k]:
            left = k
        else:
            left = np.searchsorted(t, xval) - 1

        # fill a row
        bb = _bspl.evaluate_all_bspl(t, k, xval, left)
        A[j + offset, left-k:left+1] = bb

    c = sl.solve(A, y)
    return c


def make_lsq_full_matrix(x, y, t, k=3):
    """Make the least-square spline, full matrices."""
    x, y, t = map(np.asarray, (x, y, t))
    m = x.size
    n = t.size - k - 1

    A = np.zeros((m, n), dtype=np.float_)

    for j in range(m):
        xval = x[j]
        # find interval
        if xval == t[k]:
            left = k
        else:
            left = np.searchsorted(t, xval) - 1

        # fill a row
        bb = _bspl.evaluate_all_bspl(t, k, xval, left)
        A[j, left-k:left+1] = bb

    # have observation matrix, can solve the LSQ problem
    B = np.dot(A.T, A)
    Y = np.dot(A.T, y)
    c = sl.solve(B, Y)

    return c, (A, Y)


class TestLSQ(object):
    #
    # Test make_lsq_spline
    #
    np.random.seed(1234)
    n, k = 13, 3
    x = np.sort(np.random.random(n))
    y = np.random.random(n)
    t = _augknt(np.linspace(x[0], x[-1], 7), k)

    def test_lstsq(self):
        # check LSQ construction vs a full matrix version
        x, y, t, k = self.x, self.y, self.t, self.k

        c0, AY = make_lsq_full_matrix(x, y, t, k)
        b = make_lsq_spline(x, y, t, k)

        assert_allclose(b.c, c0)
        assert_equal(b.c.shape, (t.size - k - 1,))

        # also check against numpy.lstsq
        aa, yy = AY
        c1, _, _, _ = np.linalg.lstsq(aa, y, rcond=-1)
        assert_allclose(b.c, c1)

    def test_weights(self):
        # weights = 1 is same as None
        x, y, t, k = self.x, self.y, self.t, self.k
        w = np.ones_like(x)

        b = make_lsq_spline(x, y, t, k)
        b_w = make_lsq_spline(x, y, t, k, w=w)

        assert_allclose(b.t, b_w.t, atol=1e-14)
        assert_allclose(b.c, b_w.c, atol=1e-14)
        assert_equal(b.k, b_w.k)

    def test_multiple_rhs(self):
        x, t, k, n = self.x, self.t, self.k, self.n
        y = np.random.random(size=(n, 5, 6, 7))

        b = make_lsq_spline(x, y, t, k)
        assert_equal(b.c.shape, (t.size-k-1, 5, 6, 7))

    def test_complex(self):
        # cmplx-valued `y`
        x, t, k = self.x, self.t, self.k
        yc = self.y * (1. + 2.j)

        b = make_lsq_spline(x, yc, t, k)
        b_re = make_lsq_spline(x, yc.real, t, k)
        b_im = make_lsq_spline(x, yc.imag, t, k)

        assert_allclose(b(x), b_re(x) + 1.j*b_im(x), atol=1e-15, rtol=1e-15)

    def test_int_xy(self):
        x = np.arange(10).astype(np.int_)
        y = np.arange(10).astype(np.int_)
        t = _augknt(x, k=1)
        # cython chokes on "buffer type mismatch"
        make_lsq_spline(x, y, t, k=1)

    def test_sliced_input(self):
        # cython code chokes on non C contiguous arrays
        xx = np.linspace(-1, 1, 100)

        x = xx[::3]
        y = xx[::3]
        t = _augknt(x, 1)
        make_lsq_spline(x, y, t, k=1)

    def test_checkfinite(self):
        # check_finite defaults to True; nans and such trigger a ValueError
        x = np.arange(12).astype(float)
        y = x**2
        t = _augknt(x, 3)

        for z in [np.nan, np.inf, -np.inf]:
            y[-1] = z
            assert_raises(ValueError, make_lsq_spline, x, y, t)

