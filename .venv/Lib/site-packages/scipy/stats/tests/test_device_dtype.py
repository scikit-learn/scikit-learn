import pytest
import numpy as np
from scipy import stats

from scipy._lib._array_api import xp_device, is_array_api_strict, is_torch
from scipy.stats._stats_py import _xp_mean, _xp_var

skip_xp_backends = pytest.mark.skip_xp_backends

dtypes = ['float32', 'float64']


if np.__version__ < "2":  # need NEP 50 dtype behavior
    pytest.skip(allow_module_level=True)


def get_arrays(n_arrays, *, dtype=np.float64, xp=np, shape=(30,), device=None,
               seed=84912165484321):
    rng = np.random.default_rng(seed)

    datas = []
    for i in range(n_arrays):
        data = 10*rng.random(size=shape)
        if xp.isdtype(dtype, 'complex floating'):
            data = data * 10j*rng.standard_normal(size=shape)
        data = xp.asarray(data, dtype=dtype, device=device)
        datas.append(data)

    return datas


@pytest.mark.parametrize('fun, kwargs', [(stats.gmean, {}),
                                         (stats.hmean, {}),
                                         (stats.pmean, {'p': 2}),
                                         (_xp_mean, {}),
                                         ])
@pytest.mark.parametrize('dtype', dtypes)
def test_xmean(fun, kwargs, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        x, weights = get_arrays(2, device=device, dtype=dtype, xp=xp)
        res = fun(x, weights=weights, **kwargs)
        assert xp_device(res) == xp_device(x)
        assert x.dtype == dtype


@skip_xp_backends('array_api_strict',
                  reason="special functions don't work with 'device1'")
@pytest.mark.parametrize('nargs', [1, 2])
@pytest.mark.parametrize('dtype', dtypes)
def test_entropy(nargs, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        args = get_arrays(nargs, device=device, dtype=dtype, xp=xp)
        res = stats.entropy(*args)
        assert xp_device(res) == xp_device(args[0])
        assert res.dtype == dtype


@pytest.mark.parametrize('dtype', dtypes)
def test_directional_stats(dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        x = get_arrays(1, shape=(30, 3), device=device, dtype=dtype, xp=xp)[0]
        res = stats.directional_stats(x)
        assert xp_device(res.mean_direction) == xp_device(x)
        assert res.mean_direction.dtype == dtype
        assert xp_device(res.mean_resultant_length) == xp_device(x)
        assert res.mean_resultant_length.dtype == dtype


@pytest.mark.parametrize('method', [
    "vasicek",
    "van es",
    pytest.param(
        "correa",
        marks=[
            skip_xp_backends("array_api_strict", reason="Invalid fancy indexing"),
            skip_xp_backends("dask.array", reason="Invalid fancy indexing"),
        ],
    ),
    "ebrahimi",
    "auto",
])
@pytest.mark.parametrize('dtype', dtypes)
def test_differential_entropy(method, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        values = get_arrays(1, device=device, dtype=dtype, xp=xp)[0]
        res = stats.differential_entropy(values, method=method)
        assert xp_device(res) == xp_device(values)
        assert values.dtype == dtype


@skip_xp_backends('dask.array', reason='no take_along_axis')
@pytest.mark.parametrize('method', ['inverted_cdf', 'averaged_inverted_cdf',
                                    'closest_observation', 'interpolated_inverted_cdf',
                                    'hazen', 'weibull', 'linear', 'median_unbiased',
                                    'normal_unbiased', 'harrell-davis'])
@pytest.mark.parametrize('dtype', dtypes)
def test_quantile(method, dtype, xp, devices):
    if (is_array_api_strict(xp) or is_torch(xp)) and method == 'harrell-davis':
        pytest.skip("'harrell-davis' not currently supported on GPU.")

    dtype = getattr(xp, dtype)
    for device in devices:
        values = get_arrays(1, device=device, dtype=dtype, xp=xp)[0]
        res = stats.quantile(values, 0.5, method=method)
        assert xp_device(res) == xp_device(values)
        assert values.dtype == dtype


@pytest.mark.parametrize('dtype', dtypes)
def test_boxcox_llf(dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        data = get_arrays(1, device=device, dtype=dtype, xp=xp)[0]
        res = stats.boxcox_llf(1, data)
        assert xp_device(res) == xp_device(data)
        assert data.dtype == dtype


@pytest.mark.parametrize('fun, kwargs',
    [(stats.moment, {'order': 2}),
     (stats.skew, {}),
     (stats.skew, {'bias': False}),
     (stats.kurtosis, {}),
     (stats.kurtosis, {'bias': False}),
     (stats.sem, {}),
     (stats.kstat, {'n': 1}),
     (stats.kstat, {'n': 2}),
     (stats.kstat, {'n': 3}),
     (stats.kstat, {'n': 4}),
     (stats.kstatvar, {'n': 1}),
     (stats.kstatvar, {'n': 2}),
     (stats.circmean, {}),
     (stats.circvar, {}),
     (stats.circstd, {}),
     (_xp_var, {}),
     (stats.gstd, {}),
     (stats.variation, {}),
     (stats.tmean, {'limits': (0.1, 0.9)}),
     (stats.tvar, {'limits': (0.1, 0.9)}),
     (stats.tmin, {'lowerlimit': 0.1}),
     (stats.tmax, {'upperlimit': 0.9}),
     (stats.tstd, {'limits': (0.1, 0.9)}),
     (stats.tsem, {'limits': (0.1, 0.9)}),
     ])
@pytest.mark.parametrize('dtype', dtypes)
def test_one_in_one_out(fun, kwargs, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        array = get_arrays(1, device=device, dtype=dtype, xp=xp)[0]
        res = fun(array, **kwargs)
        assert xp_device(res) == xp_device(array)
        assert res.dtype == dtype


@pytest.mark.parametrize('dtype', dtypes)
def test_describe(dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        array = get_arrays(1, shape=10, dtype=dtype, device=device, xp=xp)[0]
        res = stats.describe(array, axis=-1)

        assert xp_device(res.nobs) == xp_device(array)
        assert xp_device(res.minmax[0]) == xp_device(array)
        assert xp_device(res.minmax[1]) == xp_device(array)
        assert xp_device(res.variance) == xp_device(array)
        assert xp_device(res.skewness) == xp_device(array)
        assert xp_device(res.kurtosis) == xp_device(array)

        assert res.minmax[0].dtype == dtype
        assert res.minmax[1].dtype == dtype
        assert res.variance.dtype == dtype
        assert res.skewness.dtype == dtype
        assert res.kurtosis.dtype == dtype


@pytest.mark.parametrize('fun', [stats.zscore, stats.gzscore, stats.zmap])
@pytest.mark.parametrize('dtype', dtypes)
def test_zscore(fun, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        n = 2 if (fun == stats.zmap) else 1
        arrays = get_arrays(n, device=device, dtype=dtype, xp=xp)
        res = fun(*arrays)
        assert xp_device(res) == xp_device(arrays[0])
        assert res.dtype == dtype


@skip_xp_backends('array_api_strict', reason="special/_support_alternative_backends")
@skip_xp_backends(cpu_only=True, exceptions=['cupy', 'jax.numpy'])
@pytest.mark.parametrize('f_name', ['ttest_1samp', 'ttest_rel', 'ttest_ind', 'skewtest',
                                    'kurtosistest', 'normaltest', 'jarque_bera',
                                    'bartlett', 'pearsonr', 'chisquare',
                                    'power_divergence'])
@pytest.mark.parametrize('dtype', dtypes)
def test_hypothesis_tests(f_name, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        f = getattr(stats, f_name)

        n = 2 if f_name in {'ttest_1samp', 'ttest_rel', 'ttest_ind', 'bartlett',
                            'pearsonr', 'chisquare'} else 1

        arrays = get_arrays(n, xp=xp, dtype=dtype, device=device)
        if f_name == 'ttest_1samp':
            arrays[1] = xp.mean(arrays[1])
        if f_name == 'chisquare':
            arrays[1] = xp.sum(arrays[0]) * arrays[1] / xp.sum(arrays[1])

        res = f(*arrays)
        assert xp_device(res.statistic) == xp_device(arrays[0])
        assert xp_device(res.pvalue) == xp_device(arrays[0])
        assert res.statistic.dtype == dtype
        assert res.pvalue.dtype == dtype

        if f_name in {'ttest_1samp', 'ttest_rel', 'ttest_ind', 'pearsonr'}:
            res_ci = res.confidence_interval()
            assert xp_device(res_ci.low) == xp_device(arrays[0])
            assert xp_device(res_ci.high) == xp_device(arrays[0])
            assert res_ci.low.dtype == dtype
            assert res_ci.high.dtype == dtype


@skip_xp_backends('array_api_strict',
                  reason="special functions don't work with 'device1'")
@skip_xp_backends(cpu_only=True, exceptions=['cupy', 'jax.numpy'])
@pytest.mark.parametrize('method', ['fisher', 'pearson', 'tippett', 'stouffer',
                                    'mudholkar_george'])
@pytest.mark.parametrize('dtype', dtypes)
def test_combine_pvalues(method, dtype, xp, devices):
    dtype = getattr(xp, dtype)
    for device in devices:
        pvalues = get_arrays(1, xp=xp, dtype=dtype, device=device)[0] / 10
        res = stats.combine_pvalues(pvalues, method=method)
        assert xp_device(res.statistic) == xp_device(pvalues)
        assert xp_device(res.pvalue) == xp_device(pvalues)
        assert res.statistic.dtype == dtype
        assert res.pvalue.dtype == dtype
