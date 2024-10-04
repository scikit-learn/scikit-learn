# mypy: disable-error-code="attr-defined"
import os
import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

import scipy._lib._elementwise_iterative_method as eim
from scipy import special, stats
from scipy.integrate import quad_vec
from scipy.integrate._tanhsinh import _tanhsinh, _pair_cache, _nsum
from scipy.stats._discrete_distns import _gen_harmonic_gt1

class TestTanhSinh:

    # Test problems from [1] Section 6
    def f1(self, t):
        return t * np.log(1 + t)

    f1.ref = 0.25
    f1.b = 1

    def f2(self, t):
        return t ** 2 * np.arctan(t)

    f2.ref = (np.pi - 2 + 2 * np.log(2)) / 12
    f2.b = 1

    def f3(self, t):
        return np.exp(t) * np.cos(t)

    f3.ref = (np.exp(np.pi / 2) - 1) / 2
    f3.b = np.pi / 2

    def f4(self, t):
        a = np.sqrt(2 + t ** 2)
        return np.arctan(a) / ((1 + t ** 2) * a)

    f4.ref = 5 * np.pi ** 2 / 96
    f4.b = 1

    def f5(self, t):
        return np.sqrt(t) * np.log(t)

    f5.ref = -4 / 9
    f5.b = 1

    def f6(self, t):
        return np.sqrt(1 - t ** 2)

    f6.ref = np.pi / 4
    f6.b = 1

    def f7(self, t):
        return np.sqrt(t) / np.sqrt(1 - t ** 2)

    f7.ref = 2 * np.sqrt(np.pi) * special.gamma(3 / 4) / special.gamma(1 / 4)
    f7.b = 1

    def f8(self, t):
        return np.log(t) ** 2

    f8.ref = 2
    f8.b = 1

    def f9(self, t):
        return np.log(np.cos(t))

    f9.ref = -np.pi * np.log(2) / 2
    f9.b = np.pi / 2

    def f10(self, t):
        return np.sqrt(np.tan(t))

    f10.ref = np.pi * np.sqrt(2) / 2
    f10.b = np.pi / 2

    def f11(self, t):
        return 1 / (1 + t ** 2)

    f11.ref = np.pi / 2
    f11.b = np.inf

    def f12(self, t):
        return np.exp(-t) / np.sqrt(t)

    f12.ref = np.sqrt(np.pi)
    f12.b = np.inf

    def f13(self, t):
        return np.exp(-t ** 2 / 2)

    f13.ref = np.sqrt(np.pi / 2)
    f13.b = np.inf

    def f14(self, t):
        return np.exp(-t) * np.cos(t)

    f14.ref = 0.5
    f14.b = np.inf

    def f15(self, t):
        return np.sin(t) / t

    f15.ref = np.pi / 2
    f15.b = np.inf

    def error(self, res, ref, log=False):
        err = abs(res - ref)

        if not log:
            return err

        with np.errstate(divide='ignore'):
            return np.log10(err)

    def test_input_validation(self):
        f = self.f1

        message = '`f` must be callable.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(42, 0, f.b)

        message = '...must be True or False.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, log=2)

        message = '...must be real numbers.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 1+1j, f.b)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol='ekki')
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=pytest)

        message = '...must be non-negative and finite.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol=np.inf)

        message = '...may not be positive infinity.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=np.inf, log=True)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol=np.inf, log=True)

        message = '...must be integers.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxlevel=object())
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxfun=1+1j)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minlevel="migratory coconut")

        message = '...must be non-negative.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxlevel=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxfun=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minlevel=-1)

        message = '...must be True or False.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, preserve_shape=2)

        message = '...must be callable.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, callback='elderberry')

    @pytest.mark.parametrize("limits, ref", [
        [(0, np.inf), 0.5],  # b infinite
        [(-np.inf, 0), 0.5],  # a infinite
        [(-np.inf, np.inf), 1],  # a and b infinite
        [(np.inf, -np.inf), -1],  # flipped limits
        [(1, -1), stats.norm.cdf(-1) -  stats.norm.cdf(1)],  # flipped limits
    ])
    def test_integral_transforms(self, limits, ref):
        # Check that the integral transforms are behaving for both normal and
        # log integration
        dist = stats.norm()

        res = _tanhsinh(dist.pdf, *limits)
        assert_allclose(res.integral, ref)

        logres = _tanhsinh(dist.logpdf, *limits, log=True)
        assert_allclose(np.exp(logres.integral), ref)
        # Transformation should not make the result complex unnecessarily
        assert (np.issubdtype(logres.integral.dtype, np.floating) if ref > 0
                else np.issubdtype(logres.integral.dtype, np.complexfloating))

        assert_allclose(np.exp(logres.error), res.error, atol=1e-16)

    # 15 skipped intentionally; it's very difficult numerically
    @pytest.mark.parametrize('f_number', range(1, 15))
    def test_basic(self, f_number):
        f = getattr(self, f"f{f_number}")
        rtol = 2e-8
        res = _tanhsinh(f, 0, f.b, rtol=rtol)
        assert_allclose(res.integral, f.ref, rtol=rtol)
        if f_number not in {14}:  # mildly underestimates error here
            true_error = abs(self.error(res.integral, f.ref)/res.integral)
            assert true_error < res.error

        if f_number in {7, 10, 12}:  # succeeds, but doesn't know it
            return

        assert res.success
        assert res.status == 0

    @pytest.mark.parametrize('ref', (0.5, [0.4, 0.6]))
    @pytest.mark.parametrize('case', stats._distr_params.distcont)
    def test_accuracy(self, ref, case):
        distname, params = case
        if distname in {'dgamma', 'dweibull', 'laplace', 'kstwo'}:
            # should split up interval at first-derivative discontinuity
            pytest.skip('tanh-sinh is not great for non-smooth integrands')
        if (distname in {'studentized_range', 'levy_stable'}
                and not int(os.getenv('SCIPY_XSLOW', 0))):
            pytest.skip('This case passes, but it is too slow.')
        dist = getattr(stats, distname)(*params)
        x = dist.interval(ref)
        res = _tanhsinh(dist.pdf, *x)
        assert_allclose(res.integral, ref)

    @pytest.mark.parametrize('shape', [tuple(), (12,), (3, 4), (3, 2, 2)])
    def test_vectorization(self, shape):
        # Test for correct functionality, output shapes, and dtypes for various
        # input shapes.
        rng = np.random.default_rng(82456839535679456794)
        a = rng.random(shape)
        b = rng.random(shape)
        p = rng.random(shape)
        n = np.prod(shape)

        def f(x, p):
            f.ncall += 1
            f.feval += 1 if (x.size == n or x.ndim <=1) else x.shape[-1]
            return x**p
        f.ncall = 0
        f.feval = 0

        @np.vectorize
        def _tanhsinh_single(a, b, p):
            return _tanhsinh(lambda x: x**p, a, b)

        res = _tanhsinh(f, a, b, args=(p,))
        refs = _tanhsinh_single(a, b, p).ravel()

        attrs = ['integral', 'error', 'success', 'status', 'nfev', 'maxlevel']
        for attr in attrs:
            ref_attr = [getattr(ref, attr) for ref in refs]
            res_attr = getattr(res, attr)
            assert_allclose(res_attr.ravel(), ref_attr, rtol=1e-15)
            assert_equal(res_attr.shape, shape)

        assert np.issubdtype(res.success.dtype, np.bool_)
        assert np.issubdtype(res.status.dtype, np.integer)
        assert np.issubdtype(res.nfev.dtype, np.integer)
        assert np.issubdtype(res.maxlevel.dtype, np.integer)
        assert_equal(np.max(res.nfev), f.feval)
        # maxlevel = 2 -> 3 function calls (2 initialization, 1 work)
        assert np.max(res.maxlevel) >= 2
        assert_equal(np.max(res.maxlevel), f.ncall)

    def test_flags(self):
        # Test cases that should produce different status flags; show that all
        # can be produced simultaneously.
        def f(xs, js):
            f.nit += 1
            funcs = [lambda x: np.exp(-x**2),  # converges
                     lambda x: np.exp(x),  # reaches maxiter due to order=2
                     lambda x: np.full_like(x, np.nan)[()]]  # stops due to NaN
            res = [funcs[j](x) for x, j in zip(xs, js.ravel())]
            return res
        f.nit = 0

        args = (np.arange(3, dtype=np.int64),)
        res = _tanhsinh(f, [np.inf]*3, [-np.inf]*3, maxlevel=5, args=args)
        ref_flags = np.array([0, -2, -3])
        assert_equal(res.status, ref_flags)

    def test_flags_preserve_shape(self):
        # Same test as above but using `preserve_shape` option to simplify.
        def f(x):
            return [np.exp(-x[0]**2),  # converges
                    np.exp(x[1]),  # reaches maxiter due to order=2
                    np.full_like(x[2], np.nan)[()]]  # stops due to NaN

        res = _tanhsinh(f, [np.inf]*3, [-np.inf]*3, maxlevel=5, preserve_shape=True)
        ref_flags = np.array([0, -2, -3])
        assert_equal(res.status, ref_flags)

    def test_preserve_shape(self):
        # Test `preserve_shape` option
        def f(x):
            return np.asarray([[x, np.sin(10 * x)],
                               [np.cos(30 * x), x * np.sin(100 * x)]])

        ref = quad_vec(f, 0, 1)
        res = _tanhsinh(f, 0, 1, preserve_shape=True)
        assert_allclose(res.integral, ref[0])

    def test_convergence(self):
        # demonstrate that number of accurate digits doubles each iteration
        f = self.f1
        last_logerr = 0
        for i in range(4):
            res = _tanhsinh(f, 0, f.b, minlevel=0, maxlevel=i)
            logerr = self.error(res.integral, f.ref, log=True)
            assert (logerr < last_logerr * 2 or logerr < -15.5)
            last_logerr = logerr

    def test_options_and_result_attributes(self):
        # demonstrate that options are behaving as advertised and status
        # messages are as intended
        def f(x):
            f.calls += 1
            f.feval += np.size(x)
            return self.f2(x)
        f.ref = self.f2.ref
        f.b = self.f2.b
        default_rtol = 1e-12
        default_atol = f.ref * default_rtol  # effective default absolute tol

        # Test default options
        f.feval, f.calls = 0, 0
        ref = _tanhsinh(f, 0, f.b)
        assert self.error(ref.integral, f.ref) < ref.error < default_atol
        assert ref.nfev == f.feval
        ref.calls = f.calls  # reference number of function calls
        assert ref.success
        assert ref.status == 0

        # Test `maxlevel` equal to required max level
        # We should get all the same results
        f.feval, f.calls = 0, 0
        maxlevel = ref.maxlevel
        res = _tanhsinh(f, 0, f.b, maxlevel=maxlevel)
        res.calls = f.calls
        assert res == ref

        # Now reduce the maximum level. We won't meet tolerances.
        f.feval, f.calls = 0, 0
        maxlevel -= 1
        assert maxlevel >= 2  # can't compare errors otherwise
        res = _tanhsinh(f, 0, f.b, maxlevel=maxlevel)
        assert self.error(res.integral, f.ref) < res.error > default_atol
        assert res.nfev == f.feval < ref.nfev
        assert f.calls == ref.calls - 1
        assert not res.success
        assert res.status == eim._ECONVERR

        # `maxfun` is currently not enforced

        # # Test `maxfun` equal to required number of function evaluations
        # # We should get all the same results
        # f.feval, f.calls = 0, 0
        # maxfun = ref.nfev
        # res = _tanhsinh(f, 0, f.b, maxfun = maxfun)
        # assert res == ref
        #
        # # Now reduce `maxfun`. We won't meet tolerances.
        # f.feval, f.calls = 0, 0
        # maxfun -= 1
        # res = _tanhsinh(f, 0, f.b, maxfun=maxfun)
        # assert self.error(res.integral, f.ref) < res.error > default_atol
        # assert res.nfev == f.feval < ref.nfev
        # assert f.calls == ref.calls - 1
        # assert not res.success
        # assert res.status == 2

        # Take this result to be the new reference
        ref = res
        ref.calls = f.calls

        # Test `atol`
        f.feval, f.calls = 0, 0
        # With this tolerance, we should get the exact same result as ref
        atol = np.nextafter(ref.error, np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=0, atol=atol)
        assert res.integral == ref.integral
        assert res.error == ref.error
        assert res.nfev == f.feval == ref.nfev
        assert f.calls == ref.calls
        # Except the result is considered to be successful
        assert res.success
        assert res.status == 0

        f.feval, f.calls = 0, 0
        # With a tighter tolerance, we should get a more accurate result
        atol = np.nextafter(ref.error, -np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=0, atol=atol)
        assert self.error(res.integral, f.ref) < res.error < atol
        assert res.nfev == f.feval > ref.nfev
        assert f.calls > ref.calls
        assert res.success
        assert res.status == 0

        # Test `rtol`
        f.feval, f.calls = 0, 0
        # With this tolerance, we should get the exact same result as ref
        rtol = np.nextafter(ref.error/ref.integral, np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=rtol)
        assert res.integral == ref.integral
        assert res.error == ref.error
        assert res.nfev == f.feval == ref.nfev
        assert f.calls == ref.calls
        # Except the result is considered to be successful
        assert res.success
        assert res.status == 0

        f.feval, f.calls = 0, 0
        # With a tighter tolerance, we should get a more accurate result
        rtol = np.nextafter(ref.error/ref.integral, -np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=rtol)
        assert self.error(res.integral, f.ref)/f.ref < res.error/res.integral < rtol
        assert res.nfev == f.feval > ref.nfev
        assert f.calls > ref.calls
        assert res.success
        assert res.status == 0

    @pytest.mark.parametrize('rtol', [1e-4, 1e-14])
    def test_log(self, rtol):
        # Test equivalence of log-integration and regular integration
        dist = stats.norm()

        test_tols = dict(atol=1e-18, rtol=1e-15)

        # Positive integrand (real log-integrand)
        res = _tanhsinh(dist.logpdf, -1, 2, log=True, rtol=np.log(rtol))
        ref = _tanhsinh(dist.pdf, -1, 2, rtol=rtol)
        assert_allclose(np.exp(res.integral), ref.integral, **test_tols)
        assert_allclose(np.exp(res.error), ref.error, **test_tols)
        assert res.nfev == ref.nfev

        # Real integrand (complex log-integrand)
        def f(x):
            return -dist.logpdf(x)*dist.pdf(x)

        def logf(x):
            return np.log(dist.logpdf(x) + 0j) + dist.logpdf(x) + np.pi * 1j

        res = _tanhsinh(logf, -np.inf, np.inf, log=True)
        ref = _tanhsinh(f, -np.inf, np.inf)
        # In gh-19173, we saw `invalid` warnings on one CI platform.
        # Silencing `all` because I can't reproduce locally and don't want
        # to risk the need to run CI again.
        with np.errstate(all='ignore'):
            assert_allclose(np.exp(res.integral), ref.integral, **test_tols)
            assert_allclose(np.exp(res.error), ref.error, **test_tols)
        assert res.nfev == ref.nfev

    def test_complex(self):
        # Test integration of complex integrand
        # Finite limits
        def f(x):
            return np.exp(1j * x)

        res = _tanhsinh(f, 0, np.pi/4)
        ref = np.sqrt(2)/2 + (1-np.sqrt(2)/2)*1j
        assert_allclose(res.integral, ref)

        # Infinite limits
        dist1 = stats.norm(scale=1)
        dist2 = stats.norm(scale=2)
        def f(x):
            return dist1.pdf(x) + 1j*dist2.pdf(x)

        res = _tanhsinh(f, np.inf, -np.inf)
        assert_allclose(res.integral, -(1+1j))

    @pytest.mark.parametrize("maxlevel", range(4))
    def test_minlevel(self, maxlevel):
        # Verify that minlevel does not change the values at which the
        # integrand is evaluated or the integral/error estimates, only the
        # number of function calls
        def f(x):
            f.calls += 1
            f.feval += np.size(x)
            f.x = np.concatenate((f.x, x.ravel()))
            return self.f2(x)
        f.feval, f.calls, f.x = 0, 0, np.array([])

        ref = _tanhsinh(f, 0, self.f2.b, minlevel=0, maxlevel=maxlevel)
        ref_x = np.sort(f.x)

        for minlevel in range(0, maxlevel + 1):
            f.feval, f.calls, f.x = 0, 0, np.array([])
            options = dict(minlevel=minlevel, maxlevel=maxlevel)
            res = _tanhsinh(f, 0, self.f2.b, **options)
            # Should be very close; all that has changed is the order of values
            assert_allclose(res.integral, ref.integral, rtol=4e-16)
            # Difference in absolute errors << magnitude of integral
            assert_allclose(res.error, ref.error, atol=4e-16 * ref.integral)
            assert res.nfev == f.feval == len(f.x)
            assert f.calls == maxlevel - minlevel + 1 + 1  # 1 validation call
            assert res.status == ref.status
            assert_equal(ref_x, np.sort(f.x))

    def test_improper_integrals(self):
        # Test handling of infinite limits of integration (mixed with finite limits)
        def f(x):
            x[np.isinf(x)] = np.nan
            return np.exp(-x**2)
        a = [-np.inf, 0, -np.inf, np.inf, -20, -np.inf, -20]
        b = [np.inf, np.inf, 0, -np.inf, 20, 20, np.inf]
        ref = np.sqrt(np.pi)
        res = _tanhsinh(f, a, b)
        assert_allclose(res.integral, [ref, ref/2, ref/2, -ref, ref, ref, ref])

    @pytest.mark.parametrize("limits", ((0, 3), ([-np.inf, 0], [3, 3])))
    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_dtype(self, limits, dtype):
        # Test that dtypes are preserved
        a, b = np.asarray(limits, dtype=dtype)[()]

        def f(x):
            assert x.dtype == dtype
            return np.exp(x)

        rtol = 1e-12 if dtype == np.float64 else 1e-5
        res = _tanhsinh(f, a, b, rtol=rtol)
        assert res.integral.dtype == dtype
        assert res.error.dtype == dtype
        assert np.all(res.success)
        assert_allclose(res.integral, np.exp(b)-np.exp(a), rtol=rtol)

    def test_maxiter_callback(self):
        # Test behavior of `maxiter` parameter and `callback` interface
        a, b = -np.inf, np.inf
        def f(x):
            return np.exp(-x*x)

        minlevel, maxlevel = 0, 2
        maxiter = maxlevel - minlevel + 1
        kwargs = dict(minlevel=minlevel, maxlevel=maxlevel, rtol=1e-15)
        res = _tanhsinh(f, a, b, **kwargs)
        assert not res.success
        assert res.maxlevel == maxlevel

        def callback(res):
            callback.iter += 1
            callback.res = res
            assert hasattr(res, 'integral')
            assert res.status == 1
            if callback.iter == maxiter:
                raise StopIteration
        callback.iter = -1  # callback called once before first iteration
        callback.res = None

        del kwargs['maxlevel']
        res2 = _tanhsinh(f, a, b, **kwargs, callback=callback)
        # terminating with callback is identical to terminating due to maxiter
        # (except for `status`)
        for key in res.keys():
            if key == 'status':
                assert callback.res[key] == 1
                assert res[key] == -2
                assert res2[key] == -4
            else:
                assert res2[key] == callback.res[key] == res[key]

    def test_jumpstart(self):
        # The intermediate results at each level i should be the same as the
        # final results when jumpstarting at level i; i.e. minlevel=maxlevel=i
        a, b = -np.inf, np.inf
        def f(x):
            return np.exp(-x*x)

        def callback(res):
            callback.integrals.append(res.integral)
            callback.errors.append(res.error)
        callback.integrals = []
        callback.errors = []

        maxlevel = 4
        _tanhsinh(f, a, b, minlevel=0, maxlevel=maxlevel, callback=callback)

        integrals = []
        errors = []
        for i in range(maxlevel + 1):
            res = _tanhsinh(f, a, b, minlevel=i, maxlevel=i)
            integrals.append(res.integral)
            errors.append(res.error)

        assert_allclose(callback.integrals[1:], integrals, rtol=1e-15)
        assert_allclose(callback.errors[1:], errors, rtol=1e-15, atol=1e-16)

    def test_special_cases(self):
        # Test edge cases and other special cases

        # Test that integers are not passed to `f`
        # (otherwise this would overflow)
        def f(x):
            assert np.issubdtype(x.dtype, np.floating)
            return x ** 99

        res = _tanhsinh(f, 0, 1)
        assert res.success
        assert_allclose(res.integral, 1/100)

        # Test levels 0 and 1; error is NaN
        res = _tanhsinh(f, 0, 1, maxlevel=0)
        assert res.integral > 0
        assert_equal(res.error, np.nan)
        res = _tanhsinh(f, 0, 1, maxlevel=1)
        assert res.integral > 0
        assert_equal(res.error, np.nan)

        # Tes equal left and right integration limits
        res = _tanhsinh(f, 1, 1)
        assert res.success
        assert res.maxlevel == -1
        assert_allclose(res.integral, 0)

        # Test scalar `args` (not in tuple)
        def f(x, c):
            return x**c

        res = _tanhsinh(f, 0, 1, args=99)
        assert_allclose(res.integral, 1/100)

        # Test NaNs
        a = [np.nan, 0, 0, 0]
        b = [1, np.nan, 1, 1]
        c = [1, 1, np.nan, 1]
        res = _tanhsinh(f, a, b, args=(c,))
        assert_allclose(res.integral, [np.nan, np.nan, np.nan, 0.5])
        assert_allclose(res.error[:3], np.nan)
        assert_equal(res.status, [-3, -3, -3, 0])
        assert_equal(res.success, [False, False, False, True])
        assert_equal(res.nfev[:3], 1)

        # Test complex integral followed by real integral
        # Previously, h0 was of the result dtype. If the `dtype` were complex,
        # this could lead to complex cached abscissae/weights. If these get
        # cast to real dtype for a subsequent real integral, we would get a
        # ComplexWarning. Check that this is avoided.
        _pair_cache.xjc = np.empty(0)
        _pair_cache.wj = np.empty(0)
        _pair_cache.indices = [0]
        _pair_cache.h0 = None
        res = _tanhsinh(lambda x: x*1j, 0, 1)
        assert_allclose(res.integral, 0.5*1j)
        res = _tanhsinh(lambda x: x, 0, 1)
        assert_allclose(res.integral, 0.5)

        # Test zero-size
        shape = (0, 3)
        res = _tanhsinh(lambda x: x, 0, np.zeros(shape))
        attrs = ['integral', 'error', 'success', 'status', 'nfev', 'maxlevel']
        for attr in attrs:
            assert_equal(res[attr].shape, shape)


class TestNSum:
    rng = np.random.default_rng(5895448232066142650)
    p = rng.uniform(1, 10, size=10)

    def f1(self, k):
        # Integers are never passed to `f1`; if they were, we'd get
        # integer to negative integer power error
        return k**(-2)

    f1.ref = np.pi**2/6
    f1.a = 1
    f1.b = np.inf
    f1.args = tuple()

    def f2(self, k, p):
        return 1 / k**p

    f2.ref = special.zeta(p, 1)
    f2.a = 1
    f2.b = np.inf
    f2.args = (p,)

    def f3(self, k, p):
        return 1 / k**p

    f3.a = 1
    f3.b = rng.integers(5, 15, size=(3, 1))
    f3.ref = _gen_harmonic_gt1(f3.b, p)
    f3.args = (p,)

    def test_input_validation(self):
        f = self.f1

        message = '`f` must be callable.'
        with pytest.raises(ValueError, match=message):
            _nsum(42, f.a, f.b)

        message = '...must be True or False.'
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, log=2)

        message = '...must be real numbers.'
        with pytest.raises(ValueError, match=message):
            _nsum(f, 1+1j, f.b)
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, None)
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, step=object())
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, atol='ekki')
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, rtol=pytest)

        with np.errstate(all='ignore'):
            res = _nsum(f, [np.nan, -np.inf, np.inf], 1)
            assert np.all((res.status == -1) & np.isnan(res.sum)
                          & np.isnan(res.error) & ~res.success & res.nfev == 1)
            res = _nsum(f, 10, [np.nan, 1])
            assert np.all((res.status == -1) & np.isnan(res.sum)
                          & np.isnan(res.error) & ~res.success & res.nfev == 1)
            res = _nsum(f, 1, 10, step=[np.nan, -np.inf, np.inf, -1, 0])
            assert np.all((res.status == -1) & np.isnan(res.sum)
                          & np.isnan(res.error) & ~res.success & res.nfev == 1)

        message = '...must be non-negative and finite.'
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, rtol=-1)
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, atol=np.inf)

        message = '...may not be positive infinity.'
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, rtol=np.inf, log=True)
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, atol=np.inf, log=True)

        message = '...must be a non-negative integer.'
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, maxterms=3.5)
        with pytest.raises(ValueError, match=message):
            _nsum(f, f.a, f.b, maxterms=-2)

    @pytest.mark.parametrize('f_number', range(1, 4))
    def test_basic(self, f_number):
        f = getattr(self, f"f{f_number}")
        res = _nsum(f, f.a, f.b, args=f.args)
        assert_allclose(res.sum, f.ref)
        assert_equal(res.status, 0)
        assert_equal(res.success, True)

        with np.errstate(divide='ignore'):
            logres = _nsum(lambda *args: np.log(f(*args)),
                           f.a, f.b, log=True, args=f.args)
        assert_allclose(np.exp(logres.sum), res.sum)
        assert_allclose(np.exp(logres.error), res.error)
        assert_equal(logres.status, 0)
        assert_equal(logres.success, True)

    @pytest.mark.parametrize('maxterms', [0, 1, 10, 20, 100])
    def test_integral(self, maxterms):
        # test precise behavior of integral approximation
        f = self.f1

        def logf(x):
            return -2*np.log(x)

        def F(x):
            return -1 / x

        a = np.asarray([1, 5])[:, np.newaxis]
        b = np.asarray([20, 100, np.inf])[:, np.newaxis, np.newaxis]
        step = np.asarray([0.5, 1, 2]).reshape((-1, 1, 1, 1))
        nsteps = np.floor((b - a)/step)
        b_original = b
        b = a + nsteps*step

        k = a + maxterms*step
        # partial sum
        direct = f(a + np.arange(maxterms)*step).sum(axis=-1, keepdims=True)
        integral = (F(b) - F(k))/step  # integral approximation of remainder
        low = direct + integral + f(b)  # theoretical lower bound
        high = direct + integral + f(k)  # theoretical upper bound
        ref_sum = (low + high)/2  # _nsum uses average of the two
        ref_err = (high - low)/2  # error (assuming perfect quadrature)

        # correct reference values where number of terms < maxterms
        a, b, step = np.broadcast_arrays(a, b, step)
        for i in np.ndindex(a.shape):
            ai, bi, stepi = a[i], b[i], step[i]
            if (bi - ai)/stepi + 1 <= maxterms:
                direct = f(np.arange(ai, bi+stepi, stepi)).sum()
                ref_sum[i] = direct
                ref_err[i] = direct * np.finfo(direct).eps

        rtol = 1e-12
        res = _nsum(f, a, b_original, step=step, maxterms=maxterms, rtol=rtol)
        assert_allclose(res.sum, ref_sum, rtol=10*rtol)
        assert_allclose(res.error, ref_err, rtol=100*rtol)
        assert_equal(res.status, 0)
        assert_equal(res.success, True)

        i = ((b_original - a)/step + 1 <= maxterms)
        assert_allclose(res.sum[i], ref_sum[i], rtol=1e-15)
        assert_allclose(res.error[i], ref_err[i], rtol=1e-15)

        logres = _nsum(logf, a, b_original, step=step, log=True,
                       rtol=np.log(rtol), maxterms=maxterms)
        assert_allclose(np.exp(logres.sum), res.sum)
        assert_allclose(np.exp(logres.error), res.error)
        assert_equal(logres.status, 0)
        assert_equal(logres.success, True)

    @pytest.mark.parametrize('shape', [tuple(), (12,), (3, 4), (3, 2, 2)])
    def test_vectorization(self, shape):
        # Test for correct functionality, output shapes, and dtypes for various
        # input shapes.
        rng = np.random.default_rng(82456839535679456794)
        a = rng.integers(1, 10, size=shape)
        # when the sum can be computed directly or `maxterms` is large enough
        # to meet `atol`, there are slight differences (for good reason)
        # between vectorized call and looping.
        b = np.inf
        p = rng.random(shape) + 1
        n = np.prod(shape)

        def f(x, p):
            f.feval += 1 if (x.size == n or x.ndim <= 1) else x.shape[-1]
            return 1 / x ** p

        f.feval = 0

        @np.vectorize
        def _nsum_single(a, b, p, maxterms):
            return _nsum(lambda x: 1 / x**p, a, b, maxterms=maxterms)

        res = _nsum(f, a, b, maxterms=1000, args=(p,))
        refs = _nsum_single(a, b, p, maxterms=1000).ravel()

        attrs = ['sum', 'error', 'success', 'status', 'nfev']
        for attr in attrs:
            ref_attr = [getattr(ref, attr) for ref in refs]
            res_attr = getattr(res, attr)
            assert_allclose(res_attr.ravel(), ref_attr, rtol=1e-15)
            assert_equal(res_attr.shape, shape)

        assert np.issubdtype(res.success.dtype, np.bool_)
        assert np.issubdtype(res.status.dtype, np.integer)
        assert np.issubdtype(res.nfev.dtype, np.integer)
        assert_equal(np.max(res.nfev), f.feval)

    def test_status(self):
        f = self.f2

        p = [2, 2, 0.9, 1.1]
        a = [0, 0, 1, 1]
        b = [10, np.inf, np.inf, np.inf]
        ref = special.zeta(p, 1)

        with np.errstate(divide='ignore'):  # intentionally dividing by zero
            res = _nsum(f, a, b, args=(p,))

        assert_equal(res.success, [False, False, False, True])
        assert_equal(res.status, [-3, -3, -2, 0])
        assert_allclose(res.sum[res.success], ref[res.success])

    def test_nfev(self):
        def f(x):
            f.nfev += np.size(x)
            return 1 / x**2

        f.nfev = 0
        res = _nsum(f, 1, 10)
        assert_equal(res.nfev, f.nfev)

        f.nfev = 0
        res = _nsum(f, 1, np.inf, atol=1e-6)
        assert_equal(res.nfev, f.nfev)

    def test_inclusive(self):
        # There was an edge case off-by one bug when `_direct` was called with
        # `inclusive=True`. Check that this is resolved.
        res = _nsum(lambda k: 1 / k ** 2, [1, 4], np.inf, maxterms=500, atol=0.1)
        ref = _nsum(lambda k: 1 / k ** 2, [1, 4], np.inf)
        assert np.all(res.sum > (ref.sum - res.error))
        assert np.all(res.sum < (ref.sum + res.error))

    def test_special_case(self):
        # test equal lower/upper limit
        f = self.f1
        a = b = 2
        res = _nsum(f, a, b)
        assert_equal(res.sum, f(a))

        # Test scalar `args` (not in tuple)
        res = _nsum(self.f2, 1, np.inf, args=2)
        assert_allclose(res.sum, self.f1.ref)  # f1.ref is correct w/ args=2

        # Test 0 size input
        a = np.empty((3, 1, 1))  # arbitrary broadcastable shapes
        b = np.empty((0, 1))  # could use Hypothesis
        p = np.empty(4)  # but it's overkill
        shape = np.broadcast_shapes(a.shape, b.shape, p.shape)
        res = _nsum(self.f2, a, b, args=(p,))
        assert res.sum.shape == shape
        assert res.status.shape == shape
        assert res.nfev.shape == shape

        # Test maxterms=0
        def f(x):
            with np.errstate(divide='ignore'):
                return 1 / x

        res = _nsum(f, 0, 10, maxterms=0)
        assert np.isnan(res.sum)
        assert np.isnan(res.error)
        assert res.status == -2

        res = _nsum(f, 0, 10, maxterms=1)
        assert np.isnan(res.sum)
        assert np.isnan(res.error)
        assert res.status == -3

        # Test NaNs
        # should skip both direct and integral methods if there are NaNs
        a = [np.nan, 1, 1, 1]
        b = [np.inf, np.nan, np.inf, np.inf]
        p = [2, 2, np.nan, 2]
        res = _nsum(self.f2, a, b, args=(p,))
        assert_allclose(res.sum, [np.nan, np.nan, np.nan, self.f1.ref])
        assert_allclose(res.error[:3], np.nan)
        assert_equal(res.status, [-1, -1, -3, 0])
        assert_equal(res.success, [False, False, False, True])
        # Ideally res.nfev[2] would be 1, but `tanhsinh` has some function evals
        assert_equal(res.nfev[:2], 1)

    @pytest.mark.parametrize('dtype', [np.float32, np.float64])
    def test_dtype(self, dtype):
        def f(k):
            assert k.dtype == dtype
            return 1 / k ** np.asarray(2, dtype=dtype)[()]

        a = np.asarray(1, dtype=dtype)
        b = np.asarray([10, np.inf], dtype=dtype)
        res = _nsum(f, a, b)
        assert res.sum.dtype == dtype
        assert res.error.dtype == dtype

        rtol = 1e-12 if dtype == np.float64 else 1e-6
        ref = _gen_harmonic_gt1(b, 2)
        assert_allclose(res.sum, ref, rtol=rtol)
