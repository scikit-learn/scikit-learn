import pytest
import numpy as np

from scipy import stats
from scipy.stats._quantile import _xp_searchsorted
from scipy._lib._array_api import (
    xp_default_dtype,
    is_numpy,
    is_torch,
    is_jax,
    make_xp_test_case,
    SCIPY_ARRAY_API,
    xp_size,
    xp_copy,
)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal
from scipy._lib._util import _apply_over_batch

skip_xp_backends = pytest.mark.skip_xp_backends

lazy_xp_modules = [stats]

@_apply_over_batch(('x', 1), ('p', 1))
def quantile_reference_last_axis(x, p, nan_policy, method):
    if nan_policy == 'omit':
        x = x[~np.isnan(x)]
    p_mask = np.isnan(p)
    p = p.copy()
    p[p_mask] = 0.5
    if method == 'harrell-davis':
        # hdquantiles returns masked element if length along axis is 1 (bug)
        res = (np.full_like(p, x[0]) if x.size == 1
               else stats.mstats.hdquantiles(x, p).data)
    elif method.startswith('round'):
        res = winsor_reference_1d(np.sort(x), p, method)
    else:
        res = np.quantile(x, p, method=method)

    res = np.asarray(res)
    if nan_policy == 'propagate' and np.any(np.isnan(x)):
        res[:] = np.nan

    res[p_mask] = np.nan
    return res


@np.vectorize(excluded={0, 2})  # type: ignore[call-arg]
def winsor_reference_1d(y, p, method):
    # Adapted directly from the documentation
    # Note: `y` is the sorted data array
    n = len(y)
    if method == 'round_nearest':
        j = int(np.round(p * n) if p < 0.5 else np.round(n * p - 1))
    elif method == 'round_outward':
        j = int(np.floor(p * n) if p < 0.5 else np.ceil(n * p - 1))
    elif method == 'round_inward':
        j = int(np.ceil(p * n) if p < 0.5 else np.floor(n * p - 1))
    return y[j]


def quantile_reference(x, p, *, axis, nan_policy, keepdims, method):
    x, p = np.moveaxis(x, axis, -1), np.moveaxis(np.atleast_1d(p), axis, -1)
    res = quantile_reference_last_axis(x, p, nan_policy, method)
    res = np.moveaxis(res, -1, axis)
    if not keepdims:
        res = np.squeeze(res, axis=axis)
    return res


@make_xp_test_case(stats.quantile)
class TestQuantile:

    def test_input_validation(self, xp):
        x = xp.asarray([1, 2, 3])
        p = xp.asarray(0.5)

        message = "`x` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(xp.asarray([True, False]), p)
        with pytest.raises(ValueError):
            stats.quantile(xp.asarray([1+1j, 2]), p)

        message = "`p` must have real floating dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0, 1]))

        message = "`weights` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=xp.astype(x, xp.complex64))

        message = "`axis` must be an integer or None."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=0.5)
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=(0, -1))

        message = "`axis` is not compatible with the shapes of the inputs."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=2)

        if not is_jax(xp):  # no data-dependent input validation for lazy arrays
            message = "The input contains nan values"
            with pytest.raises(ValueError, match=message):
                stats.quantile(xp.asarray([xp.nan, 1, 2]), p, nan_policy='raise')

        message = "`method` must be one of..."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, method='a duck')

        message = "`method='harrell-davis'` does not support `weights`."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=x, method='harrell-davis')

        message = "`method='round_nearest'` does not support `weights`."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=x, method='round_nearest')

        message = "If specified, `keepdims` must be True or False."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, keepdims=42)

        message = "`keepdims` may be False only if the length of `p` along `axis` is 1."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0.5, 0.6]), keepdims=False)


    def _get_weights_x_rep(self, x, axis, rng):
        x = np.swapaxes(x, axis, -1)
        ndim = x.ndim
        x = np.atleast_2d(x)
        counts = rng.integers(10, size=x.shape[-1], dtype=np.int32)
        x_rep = []
        weights = []
        for x_ in x:
            counts_ = rng.permuted(counts)
            x_rep.append(np.repeat(x_, counts_))
            weights.append(counts_)
        x_rep, weights = np.stack(x_rep), np.stack(weights)
        if ndim < 2:
            x_rep, weights = np.squeeze(x_rep, axis=0), np.squeeze(weights, axis=0)
        x_rep, weights = np.swapaxes(x_rep, -1, axis), np.swapaxes(weights, -1, axis)
        weights = np.asarray(weights, dtype=x.dtype)
        return weights, x_rep


    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy'])
    @pytest.mark.parametrize('method',
         ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
          'hazen', 'interpolated_inverted_cdf', 'linear',
          'median_unbiased', 'normal_unbiased', 'weibull',
          'harrell-davis', 'round_nearest', 'round_outward', 'round_inward',
          '_lower', '_higher', '_midpoint', '_nearest'])
    @pytest.mark.parametrize('shape_x, shape_p, axis',
        [(10, None, -1), (10, 10, -1), (10, (2, 3), -1), ((10, 2), None, 0)])
    @pytest.mark.parametrize('weights', [False, True])
    def test_against_reference(self, method, shape_x, shape_p, axis, weights, xp):
        # Test all methods with various data shapes
        if weights and (method.startswith('_') or method.startswith('round')
                        or method=='harrell-davis'):
            pytest.skip('`weights` not supported by private (legacy) methods.')
        dtype = xp_default_dtype(xp)
        rng = np.random.default_rng(23458924568734956)
        x = rng.random(size=shape_x)
        p = rng.random(size=shape_p)

        if weights:
            weights, x_rep = self._get_weights_x_rep(x, axis, rng)
        else:
            weights, x_rep = None, x

        ref = quantile_reference(
            x_rep, p, method=method[1:] if method.startswith('_') else method,
            axis=axis, nan_policy='propagate', keepdims=shape_p is not None)

        x, p = xp.asarray(x, dtype=dtype), xp.asarray(p, dtype=dtype)
        weights = weights if weights is None else xp.asarray(weights, dtype=dtype)
        res = stats.quantile(x, p, method=method, weights=weights, axis=axis)

        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.filterwarnings("ignore:torch.searchsorted:UserWarning")
    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy', 'jax.numpy'])
    @pytest.mark.parametrize('axis', [0, 1])
    @pytest.mark.parametrize('keepdims', [False, True])
    @pytest.mark.parametrize('nan_policy', ['omit', 'propagate', 'marray'])
    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('method', ['linear', 'harrell-davis', 'round_nearest'])
    @pytest.mark.parametrize('weights', [False, True])
    def test_against_reference_2(self, axis, keepdims, nan_policy,
                                 dtype, method, weights, xp):
        # Test some methods with various combinations of arguments
        if is_jax(xp) and nan_policy == 'marray':  # mdhaber/marray#146
            pytest.skip("`marray` currently incompatible with JAX")
        if weights and method in {'harrell-davis', 'round_nearest'}:
            pytest.skip("These methods don't yet support weights")
        rng = np.random.default_rng(23458924568734956)
        shape = (5, 6)
        x = rng.random(size=shape).astype(dtype)
        p = rng.random(size=shape).astype(dtype)
        mask = rng.random(size=shape) > 0.8
        assert np.any(mask)
        x[mask] = np.nan
        if not keepdims:
            p = np.mean(p, axis=axis, keepdims=True)

        # inject p = 0 and p = 1 to test edge cases
        # Currently would fail with CuPy/JAX (cupy/cupy#8934, jax-ml/jax#21900);
        # remove the `if` when those are resolved.
        if is_numpy(xp):
            p0 = p.ravel()
            p0[1] = 0.
            p0[-2] = 1.

        dtype = getattr(xp, dtype)

        if weights:
            weights, x_rep = self._get_weights_x_rep(x, axis, rng)
            weights = weights if weights is None else xp.asarray(weights)
        else:
            weights, x_rep = None, x

        if nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            if weights is not None:
                pytest.skip("MArray is not yet compatible with weights")
            marray = pytest.importorskip('marray')
            kwargs = dict(axis=axis, keepdims=keepdims, method=method)
            mxp = marray._get_namespace(xp)
            x_mp = mxp.asarray(x, mask=mask)
            weights = weights if weights is None else mxp.asarray(weights)
            res = stats.quantile(x_mp, mxp.asarray(p), weights=weights, **kwargs)
            ref = quantile_reference(x_rep, p, nan_policy='omit', **kwargs)
            xp_assert_close(res.data, xp.asarray(ref, dtype=dtype))
            return

        kwargs = dict(axis=axis, keepdims=keepdims,
                      nan_policy=nan_policy, method=method)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), weights=weights, **kwargs)
        ref = quantile_reference(x_rep, p, **kwargs)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    def test_integer_input_output_dtype(self, xp):
        res = stats.quantile(xp.arange(10, dtype=xp.int64), 0.5)
        assert res.dtype == xp_default_dtype(xp)

    @pytest.mark.parametrize('x, p, ref, kwargs',
        [([], 0.5, np.nan, {}),
         ([1, 2, 3], [-1, 0, 1, 1.5, np.nan], [np.nan, 1, 3, np.nan, np.nan], {}),
         ([1, 2, 3], [], [], {}),
         ([[np.nan, 2]], 0.5, [np.nan, 2], {'nan_policy': 'omit'}),
         ([[], []], 0.5, np.full(2, np.nan), {'axis': -1}),
         ([[], []], 0.5, np.zeros((0,)), {'axis': 0, 'keepdims': False}),
         ([[], []], 0.5, np.zeros((1, 0)), {'axis': 0, 'keepdims': True}),
         ([], [0.5, 0.6], np.full(2, np.nan), {}),
         (np.arange(1, 28).reshape((3, 3, 3)), 0.5, [[[14.]]],
          {'axis': None, 'keepdims': True}),
         ([[1, 2], [3, 4]], [0.25, 0.5, 0.75], [[1.75, 2.5, 3.25]],
          {'axis': None, 'keepdims': True}),
         # Known issue:
         # ([1, 2, 3], 0.5, 2., {'weights': [0, 0, 0]})
         # See https://github.com/scipy/scipy/pull/23941#issuecomment-3503554361
         ])
    def test_edge_cases(self, x, p, ref, kwargs, xp):
        default_dtype = xp_default_dtype(xp)
        x, p, ref = xp.asarray(x), xp.asarray(p), xp.asarray(ref, dtype=default_dtype)
        res = stats.quantile(x, p, **kwargs)
        xp_assert_equal(res, ref)

    @pytest.mark.parametrize('axis', [0, 1, 2])
    @pytest.mark.parametrize('keepdims', [False, True])
    def test_size_0(self, axis, keepdims, xp):
        shape = [3, 4, 0]
        out_shape = shape.copy()
        if keepdims:
            out_shape[axis] = 1
        else:
            out_shape.pop(axis)
        res = stats.quantile(xp.zeros(tuple(shape)), 0.5, axis=axis, keepdims=keepdims)
        assert res.shape == tuple(out_shape)

    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
         '_lower', '_higher', '_midpoint', '_nearest'])
    def test_transition(self, method, xp):
        # test that values of discontinuous estimators are correct when
        # p*n + m - 1 is integral.
        if method == 'closest_observation' and np.__version__ < '2.0.1':
            pytest.skip('Bug in np.quantile (numpy/numpy#26656) fixed in 2.0.1')
        x = np.arange(8., dtype=np.float64)
        p = np.arange(0, 1.03125, 0.03125)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), method=method)
        ref = np.quantile(x, p, method=method[1:] if method.startswith('_') else method)
        xp_assert_equal(res, xp.asarray(ref, dtype=xp.float64))

    @pytest.mark.parametrize('zero_weights', [False, True])
    def test_weights_against_numpy(self, zero_weights, xp):
        if is_numpy(xp) and xp.__version__ < "2.0":
            pytest.skip('`weights` not supported by NumPy < 2.0.')
        dtype = xp_default_dtype(xp)
        rng = np.random.default_rng(85468924398205602)
        method = 'inverted_cdf'
        x = rng.random(size=100)
        weights = rng.random(size=100)
        if zero_weights:
            weights[weights < 0.5] = 0
        p = np.linspace(0., 1., 300)
        res = stats.quantile(xp.asarray(x, dtype=dtype), xp.asarray(p, dtype=dtype),
                             method=method, weights=xp.asarray(weights, dtype=dtype))
        ref = np.quantile(x, p, method=method, weights=weights)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation', 'hazen',
         'interpolated_inverted_cdf', 'linear','median_unbiased', 'normal_unbiased',
         'weibull'])
    def test_zero_weights(self, method, xp):
        rng = np.random.default_rng(85468924398205602)

        # test 1-D versus eliminating zero-weighted values
        n = 100
        x = xp.asarray(rng.random(size=n))
        x0 = xp_copy(x)
        p = xp.asarray(rng.random(size=n))
        i_zero = xp.asarray(rng.random(size=n) < 0.1)
        weights = xp.asarray(rng.random(size=n))
        weights = xp.where(i_zero, 0., weights)
        res = stats.quantile(x, p, weights=weights, method=method)
        ref = stats.quantile(x[~i_zero], p, weights=weights[~i_zero], method=method)
        xp_assert_close(res, ref)
        xp_assert_equal(x, x0)  # no input mutation

        # test multi-D versus `nan_policy='omit'`
        shape = (5, 100)
        x = xp.asarray(rng.random(size=shape))
        x0 = xp_copy(x)
        p = xp.asarray(rng.random(size=shape))
        i_zero = xp.asarray(rng.random(size=shape) < 0.1)
        weights = xp.asarray(rng.random(size=shape))
        x_nanned = xp.where(i_zero, xp.nan, x)
        weights_zeroed = xp.where(i_zero, 0., weights)
        res = stats.quantile(x, p, weights=weights_zeroed, method=method, axis=-1)
        ref = stats.quantile(x_nanned, p, weights=weights,
                             nan_policy='omit', method=method, axis=-1)
        xp_assert_close(res, ref)
        xp_assert_equal(x, x0)  # no input mutation

    @pytest.mark.filterwarnings("ignore:torch.searchsorted:UserWarning")
    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation', 'hazen',
         'interpolated_inverted_cdf', 'linear','median_unbiased', 'normal_unbiased',
         'weibull'])
    @pytest.mark.parametrize('shape', [50, (50, 3)])
    def test_unity_weights(self, method, shape, xp):
        # Check that result is unchanged if all weights are `1.0`
        rng = np.random.default_rng(28546892439820560)
        x = xp.asarray(rng.random(size=shape))
        p = xp.asarray(rng.random(size=shape))
        weights = xp.ones_like(x)
        res = stats.quantile(x, p, weights=weights, method=method)
        ref = stats.quantile(x, p, method=method)
        xp_assert_close(res, ref)


@_apply_over_batch(('a', 1), ('v', 1))
def np_searchsorted(a, v, side):
    return np.searchsorted(a, v, side=side)


@make_xp_test_case(_xp_searchsorted)
class Test_XPSearchsorted:
    @pytest.mark.parametrize('side', ['left', 'right'])
    @pytest.mark.parametrize('ties', [False, True])
    @pytest.mark.parametrize('shape', [0, 1, 2, 10, 11, 1000, 10001,
                                       (2, 0), (0, 2), (2, 10), (2, 3, 11)])
    @pytest.mark.parametrize('nans_x', [False, True])
    @pytest.mark.parametrize('infs_x', [False, True])
    def test_nd(self, side, ties, shape, nans_x, infs_x, xp):
        if nans_x and is_torch(xp):
            pytest.skip('torch sorts NaNs differently')
        rng = np.random.default_rng(945298725498274853)
        if ties:
            x = rng.integers(5, size=shape)
        else:
            x = rng.random(shape)
        # float32 is to accommodate JAX - nextafter with `float64` is too small?
        x = np.asarray(x, dtype=np.float32)
        xr = np.nextafter(x, np.inf)
        xl = np.nextafter(x, -np.inf)
        x_ = np.asarray([-np.inf, np.inf, np.nan])
        x_ = np.broadcast_to(x_, x.shape[:-1] + (3,))
        y = rng.permuted(np.concatenate((xl, x, xr, x_), axis=-1), axis=-1)
        if nans_x:
            mask = rng.random(shape) < 0.1
            x[mask] = np.nan
        if infs_x:
            mask = rng.random(shape) < 0.1
            x[mask] = -np.inf
            mask = rng.random(shape) > 0.9
            x[mask] = np.inf
        x = np.sort(x, axis=-1)
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        xp_default_int = xp.asarray(1).dtype
        if xp_size(x) == 0 and x.ndim > 0 and x.shape[-1] != 0:
            ref = xp.empty(x.shape[:-1] + (y.shape[-1],), dtype=xp_default_int)
        else:
            ref = xp.asarray(np_searchsorted(x, y, side=side), dtype=xp_default_int)
        x, y = xp.asarray(x), xp.asarray(y)
        res = _xp_searchsorted(x, y, side=side)
        xp_assert_equal(res, ref)
