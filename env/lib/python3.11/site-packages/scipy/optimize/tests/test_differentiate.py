import pytest

import numpy as np
from numpy.testing import assert_array_less, assert_allclose, assert_equal

import scipy._lib._elementwise_iterative_method as eim
from scipy import stats, optimize
from scipy.optimize._differentiate import (_differentiate as differentiate,
                                           _jacobian as jacobian, _EERRORINCREASE)

class TestDifferentiate:

    def f(self, x):
        return stats.norm().cdf(x)

    @pytest.mark.parametrize('x', [0.6, np.linspace(-0.05, 1.05, 10)])
    def test_basic(self, x):
        # Invert distribution CDF and compare against distribution `ppf`
        res = differentiate(self.f, x)
        ref = stats.norm().pdf(x)
        np.testing.assert_allclose(res.df, ref)
        # This would be nice, but doesn't always work out. `error` is an
        # estimate, not a bound.
        assert_array_less(abs(res.df - ref), res.error)
        assert res.x.shape == ref.shape

    @pytest.mark.parametrize('case', stats._distr_params.distcont)
    def test_accuracy(self, case):
        distname, params = case
        dist = getattr(stats, distname)(*params)
        x = dist.median() + 0.1
        res = differentiate(dist.cdf, x)
        ref = dist.pdf(x)
        assert_allclose(res.df, ref, atol=1e-10)

    @pytest.mark.parametrize('order', [1, 6])
    @pytest.mark.parametrize('shape', [tuple(), (12,), (3, 4), (3, 2, 2)])
    def test_vectorization(self, order, shape):
        # Test for correct functionality, output shapes, and dtypes for various
        # input shapes.
        x = np.linspace(-0.05, 1.05, 12).reshape(shape) if shape else 0.6
        n = np.size(x)

        @np.vectorize
        def _differentiate_single(x):
            return differentiate(self.f, x, order=order)

        def f(x, *args, **kwargs):
            f.nit += 1
            f.feval += 1 if (x.size == n or x.ndim <=1) else x.shape[-1]
            return self.f(x, *args, **kwargs)
        f.nit = -1
        f.feval = 0

        res = differentiate(f, x, order=order)
        refs = _differentiate_single(x).ravel()

        ref_x = [ref.x for ref in refs]
        assert_allclose(res.x.ravel(), ref_x)
        assert_equal(res.x.shape, shape)

        ref_df = [ref.df for ref in refs]
        assert_allclose(res.df.ravel(), ref_df)
        assert_equal(res.df.shape, shape)

        ref_error = [ref.error for ref in refs]
        assert_allclose(res.error.ravel(), ref_error, atol=1e-12)
        assert_equal(res.error.shape, shape)

        ref_success = [ref.success for ref in refs]
        assert_equal(res.success.ravel(), ref_success)
        assert_equal(res.success.shape, shape)
        assert np.issubdtype(res.success.dtype, np.bool_)

        ref_flag = [ref.status for ref in refs]
        assert_equal(res.status.ravel(), ref_flag)
        assert_equal(res.status.shape, shape)
        assert np.issubdtype(res.status.dtype, np.integer)

        ref_nfev = [ref.nfev for ref in refs]
        assert_equal(res.nfev.ravel(), ref_nfev)
        assert_equal(np.max(res.nfev), f.feval)
        assert_equal(res.nfev.shape, res.x.shape)
        assert np.issubdtype(res.nfev.dtype, np.integer)

        ref_nit = [ref.nit for ref in refs]
        assert_equal(res.nit.ravel(), ref_nit)
        assert_equal(np.max(res.nit), f.nit)
        assert_equal(res.nit.shape, res.x.shape)
        assert np.issubdtype(res.nit.dtype, np.integer)

    def test_flags(self):
        # Test cases that should produce different status flags; show that all
        # can be produced simultaneously.
        rng = np.random.default_rng(5651219684984213)
        def f(xs, js):
            f.nit += 1
            funcs = [lambda x: x - 2.5,  # converges
                     lambda x: np.exp(x)*rng.random(),  # error increases
                     lambda x: np.exp(x),  # reaches maxiter due to order=2
                     lambda x: np.full_like(x, np.nan)[()]]  # stops due to NaN
            res = [funcs[j](x) for x, j in zip(xs, js.ravel())]
            return res
        f.nit = 0

        args = (np.arange(4, dtype=np.int64),)
        res = differentiate(f, [1]*4, rtol=1e-14, order=2, args=args)

        ref_flags = np.array([eim._ECONVERGED,
                              _EERRORINCREASE,
                              eim._ECONVERR,
                              eim._EVALUEERR])
        assert_equal(res.status, ref_flags)

    def test_flags_preserve_shape(self):
        # Same test as above but using `preserve_shape` option to simplify.
        rng = np.random.default_rng(5651219684984213)
        def f(x):
            return [x - 2.5,  # converges
                    np.exp(x)*rng.random(),  # error increases
                    np.exp(x),  # reaches maxiter due to order=2
                    np.full_like(x, np.nan)[()]]  # stops due to NaN

        res = differentiate(f, 1, rtol=1e-14, order=2, preserve_shape=True)

        ref_flags = np.array([eim._ECONVERGED,
                              _EERRORINCREASE,
                              eim._ECONVERR,
                              eim._EVALUEERR])
        assert_equal(res.status, ref_flags)

    def test_preserve_shape(self):
        # Test `preserve_shape` option
        def f(x):
            return [x, np.sin(3*x), x+np.sin(10*x), np.sin(20*x)*(x-1)**2]

        x = 0
        ref = [1, 3*np.cos(3*x), 1+10*np.cos(10*x),
               20*np.cos(20*x)*(x-1)**2 + 2*np.sin(20*x)*(x-1)]
        res = differentiate(f, x, preserve_shape=True)
        assert_allclose(res.df, ref)

    def test_convergence(self):
        # Test that the convergence tolerances behave as expected
        dist = stats.norm()
        x = 1
        f = dist.cdf
        ref = dist.pdf(x)
        kwargs0 = dict(atol=0, rtol=0, order=4)

        kwargs = kwargs0.copy()
        kwargs['atol'] = 1e-3
        res1 = differentiate(f, x, **kwargs)
        assert_array_less(abs(res1.df - ref), 1e-3)
        kwargs['atol'] = 1e-6
        res2 = differentiate(f, x, **kwargs)
        assert_array_less(abs(res2.df - ref), 1e-6)
        assert_array_less(abs(res2.df - ref), abs(res1.df - ref))

        kwargs = kwargs0.copy()
        kwargs['rtol'] = 1e-3
        res1 = differentiate(f, x, **kwargs)
        assert_array_less(abs(res1.df - ref), 1e-3 * np.abs(ref))
        kwargs['rtol'] = 1e-6
        res2 = differentiate(f, x, **kwargs)
        assert_array_less(abs(res2.df - ref), 1e-6 * np.abs(ref))
        assert_array_less(abs(res2.df - ref), abs(res1.df - ref))

    def test_step_parameters(self):
        # Test that step factors have the expected effect on accuracy
        dist = stats.norm()
        x = 1
        f = dist.cdf
        ref = dist.pdf(x)

        res1 = differentiate(f, x, initial_step=0.5, maxiter=1)
        res2 = differentiate(f, x, initial_step=0.05, maxiter=1)
        assert abs(res2.df - ref) < abs(res1.df - ref)

        res1 = differentiate(f, x, step_factor=2, maxiter=1)
        res2 = differentiate(f, x, step_factor=20, maxiter=1)
        assert abs(res2.df - ref) < abs(res1.df - ref)

        # `step_factor` can be less than 1: `initial_step` is the minimum step
        kwargs = dict(order=4, maxiter=1, step_direction=0)
        res = differentiate(f, x, initial_step=0.5, step_factor=0.5, **kwargs)
        ref = differentiate(f, x, initial_step=1, step_factor=2, **kwargs)
        assert_allclose(res.df, ref.df, rtol=5e-15)

        # This is a similar test for one-sided difference
        kwargs = dict(order=2, maxiter=1, step_direction=1)
        res = differentiate(f, x, initial_step=1, step_factor=2, **kwargs)
        ref = differentiate(f, x, initial_step=1/np.sqrt(2), step_factor=0.5,
                                   **kwargs)
        assert_allclose(res.df, ref.df, rtol=5e-15)

        kwargs['step_direction'] = -1
        res = differentiate(f, x, initial_step=1, step_factor=2, **kwargs)
        ref = differentiate(f, x, initial_step=1/np.sqrt(2), step_factor=0.5,
                                   **kwargs)
        assert_allclose(res.df, ref.df, rtol=5e-15)

    def test_step_direction(self):
        # test that `step_direction` works as expected
        def f(x):
            y = np.exp(x)
            y[(x < 0) + (x > 2)] = np.nan
            return y

        x = np.linspace(0, 2, 10)
        step_direction = np.zeros_like(x)
        step_direction[x < 0.6], step_direction[x > 1.4] = 1, -1
        res = differentiate(f, x, step_direction=step_direction)
        assert_allclose(res.df, np.exp(x))
        assert np.all(res.success)

    def test_vectorized_step_direction_args(self):
        # test that `step_direction` and `args` are vectorized properly
        def f(x, p):
            return x ** p

        def df(x, p):
            return p * x ** (p - 1)

        x = np.array([1, 2, 3, 4]).reshape(-1, 1, 1)
        hdir = np.array([-1, 0, 1]).reshape(1, -1, 1)
        p = np.array([2, 3]).reshape(1, 1, -1)
        res = differentiate(f, x, step_direction=hdir, args=(p,))
        ref = np.broadcast_to(df(x, p), res.df.shape)
        assert_allclose(res.df, ref)

    def test_maxiter_callback(self):
        # Test behavior of `maxiter` parameter and `callback` interface
        x = 0.612814
        dist = stats.norm()
        maxiter = 3

        def f(x):
            res = dist.cdf(x)
            return res

        default_order = 8
        res = differentiate(f, x, maxiter=maxiter, rtol=1e-15)
        assert not np.any(res.success)
        assert np.all(res.nfev == default_order + 1 + (maxiter - 1)*2)
        assert np.all(res.nit == maxiter)

        def callback(res):
            callback.iter += 1
            callback.res = res
            assert hasattr(res, 'x')
            assert res.df not in callback.dfs
            callback.dfs.add(res.df)
            assert res.status == eim._EINPROGRESS
            if callback.iter == maxiter:
                raise StopIteration
        callback.iter = -1  # callback called once before first iteration
        callback.res = None
        callback.dfs = set()

        res2 = differentiate(f, x, callback=callback, rtol=1e-15)
        # terminating with callback is identical to terminating due to maxiter
        # (except for `status`)
        for key in res.keys():
            if key == 'status':
                assert res[key] == eim._ECONVERR
                assert callback.res[key] == eim._EINPROGRESS
                assert res2[key] == eim._ECALLBACK
            else:
                assert res2[key] == callback.res[key] == res[key]

    @pytest.mark.parametrize("hdir", (-1, 0, 1))
    @pytest.mark.parametrize("x", (0.65, [0.65, 0.7]))
    @pytest.mark.parametrize("dtype", (np.float16, np.float32, np.float64))
    def test_dtype(self, hdir, x, dtype):
        # Test that dtypes are preserved
        x = np.asarray(x, dtype=dtype)[()]

        def f(x):
            assert x.dtype == dtype
            return np.exp(x)

        def callback(res):
            assert res.x.dtype == dtype
            assert res.df.dtype == dtype
            assert res.error.dtype == dtype

        res = differentiate(f, x, order=4, step_direction=hdir,
                                   callback=callback)
        assert res.x.dtype == dtype
        assert res.df.dtype == dtype
        assert res.error.dtype == dtype
        eps = np.finfo(dtype).eps
        assert_allclose(res.df, np.exp(res.x), rtol=np.sqrt(eps))

    def test_input_validation(self):
        # Test input validation for appropriate error messages

        message = '`func` must be callable.'
        with pytest.raises(ValueError, match=message):
            differentiate(None, 1)

        message = 'Abscissae and function output must be real numbers.'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, -4+1j)

        message = "When `preserve_shape=False`, the shape of the array..."
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: [1, 2, 3], [-2, -3])

        message = 'Tolerances and step parameters must be non-negative...'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, atol=-1)
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, rtol='ekki')
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, initial_step=None)
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, step_factor=object())

        message = '`maxiter` must be a positive integer.'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, maxiter=1.5)
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, maxiter=0)

        message = '`order` must be a positive integer'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, order=1.5)
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, order=0)

        message = '`preserve_shape` must be True or False.'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, preserve_shape='herring')

        message = '`callback` must be callable.'
        with pytest.raises(ValueError, match=message):
            differentiate(lambda x: x, 1, callback='shrubbery')

    def test_special_cases(self):
        # Test edge cases and other special cases

        # Test that integers are not passed to `f`
        # (otherwise this would overflow)
        def f(x):
            assert np.issubdtype(x.dtype, np.floating)
            return x ** 99 - 1

        res = differentiate(f, 7, rtol=1e-10)
        assert res.success
        assert_allclose(res.df, 99*7.**98)

        # Test that if success is achieved in the correct number
        # of iterations if function is a polynomial. Ideally, all polynomials
        # of order 0-2 would get exact result with 0 refinement iterations,
        # all polynomials of order 3-4 would be differentiated exactly after
        # 1 iteration, etc. However, it seems that _differentiate needs an
        # extra iteration to detect convergence based on the error estimate.

        for n in range(6):
            x = 1.5
            def f(x):
                return 2*x**n

            ref = 2*n*x**(n-1)

            res = differentiate(f, x, maxiter=1, order=max(1, n))
            assert_allclose(res.df, ref, rtol=1e-15)
            assert_equal(res.error, np.nan)

            res = differentiate(f, x, order=max(1, n))
            assert res.success
            assert res.nit == 2
            assert_allclose(res.df, ref, rtol=1e-15)

        # Test scalar `args` (not in tuple)
        def f(x, c):
            return c*x - 1

        res = differentiate(f, 2, args=3)
        assert_allclose(res.df, 3)

    @pytest.mark.xfail
    @pytest.mark.parametrize("case", (  # function, evaluation point
        (lambda x: (x - 1) ** 3, 1),
        (lambda x: np.where(x > 1, (x - 1) ** 5, (x - 1) ** 3), 1)
    ))
    def test_saddle_gh18811(self, case):
        # With default settings, _differentiate will not always converge when
        # the true derivative is exactly zero. This tests that specifying a
        # (tight) `atol` alleviates the problem. See discussion in gh-18811.
        atol = 1e-16
        res = differentiate(*case, step_direction=[-1, 0, 1], atol=atol)
        assert np.all(res.success)
        assert_allclose(res.df, 0, atol=atol)


class TestJacobian:

    # Example functions and Jacobians from Wikipedia:
    # https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant#Examples

    def f1(z):
        x, y = z
        return [x ** 2 * y, 5 * x + np.sin(y)]

    def df1(z):
        x, y = z
        return [[2 * x * y, x ** 2], [np.full_like(x, 5), np.cos(y)]]

    f1.mn = 2, 2  # type: ignore[attr-defined]
    f1.ref = df1  # type: ignore[attr-defined]

    def f2(z):
        r, phi = z
        return [r * np.cos(phi), r * np.sin(phi)]

    def df2(z):
        r, phi = z
        return [[np.cos(phi), -r * np.sin(phi)],
                [np.sin(phi), r * np.cos(phi)]]

    f2.mn = 2, 2  # type: ignore[attr-defined]
    f2.ref = df2  # type: ignore[attr-defined]

    def f3(z):
        r, phi, th = z
        return [r * np.sin(phi) * np.cos(th), r * np.sin(phi) * np.sin(th),
                r * np.cos(phi)]

    def df3(z):
        r, phi, th = z
        return [[np.sin(phi) * np.cos(th), r * np.cos(phi) * np.cos(th),
                 -r * np.sin(phi) * np.sin(th)],
                [np.sin(phi) * np.sin(th), r * np.cos(phi) * np.sin(th),
                 r * np.sin(phi) * np.cos(th)],
                [np.cos(phi), -r * np.sin(phi), np.zeros_like(r)]]

    f3.mn = 3, 3  # type: ignore[attr-defined]
    f3.ref = df3  # type: ignore[attr-defined]

    def f4(x):
        x1, x2, x3 = x
        return [x1, 5 * x3, 4 * x2 ** 2 - 2 * x3, x3 * np.sin(x1)]

    def df4(x):
        x1, x2, x3 = x
        one = np.ones_like(x1)
        return [[one, 0 * one, 0 * one],
                [0 * one, 0 * one, 5 * one],
                [0 * one, 8 * x2, -2 * one],
                [x3 * np.cos(x1), 0 * one, np.sin(x1)]]

    f4.mn = 3, 4  # type: ignore[attr-defined]
    f4.ref = df4  # type: ignore[attr-defined]

    def f5(x):
        x1, x2, x3 = x
        return [5 * x2, 4 * x1 ** 2 - 2 * np.sin(x2 * x3), x2 * x3]

    def df5(x):
        x1, x2, x3 = x
        one = np.ones_like(x1)
        return [[0 * one, 5 * one, 0 * one],
                [8 * x1, -2 * x3 * np.cos(x2 * x3), -2 * x2 * np.cos(x2 * x3)],
                [0 * one, x3, x2]]

    f5.mn = 3, 3  # type: ignore[attr-defined]
    f5.ref = df5  # type: ignore[attr-defined]

    rosen = optimize.rosen
    rosen.mn = 5, 1  # type: ignore[attr-defined]
    rosen.ref = optimize.rosen_der  # type: ignore[attr-defined]

    @pytest.mark.parametrize('size', [(), (6,), (2, 3)])
    @pytest.mark.parametrize('func', [f1, f2, f3, f4, f5, rosen])
    def test_examples(self, size, func):
        rng = np.random.default_rng(458912319542)
        m, n = func.mn
        x = rng.random(size=(m,) + size)
        res = jacobian(func, x).df
        ref = func.ref(x)
        np.testing.assert_allclose(res, ref, atol=1e-10)

    def test_iv(self):
        # Test input validation
        message = "Argument `x` must be at least 1-D."
        with pytest.raises(ValueError, match=message):
            jacobian(np.sin, 1, atol=-1)

        # Confirm that other parameters are being passed to `_derivative`,
        # which raises an appropriate error message.
        x = np.ones(3)
        func = optimize.rosen
        message = 'Tolerances and step parameters must be non-negative scalars.'
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, atol=-1)
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, rtol=-1)
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, initial_step=-1)
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, step_factor=-1)

        message = '`order` must be a positive integer.'
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, order=-1)

        message = '`maxiter` must be a positive integer.'
        with pytest.raises(ValueError, match=message):
            jacobian(func, x, maxiter=-1)
