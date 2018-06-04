import pytest

scipy = pytest.importorskip('scipy')
import numpy as np
import dask.array as da
from dask.array.utils import assert_eq
from dask.delayed import Delayed
import dask.array.stats
from dask.array.utils import allclose


@pytest.mark.parametrize('kind, kwargs', [
    ('skew', {}),
    ('kurtosis', {}),
    ('kurtosis', {'fisher': False}),
])
def test_measures(kind, kwargs):
    x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)
    dfunc = getattr(dask.array.stats, kind)
    sfunc = getattr(scipy.stats, kind)

    expected = sfunc(x, **kwargs)
    result = dfunc(y, **kwargs)
    assert_eq(result, expected)
    assert isinstance(result, da.Array)


def test_bias_raises():
    x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)

    with pytest.raises(NotImplementedError):
        dask.array.stats.skew(y, bias=False)

    with pytest.raises(NotImplementedError):
        dask.array.stats.kurtosis(y, bias=False)


@pytest.mark.parametrize('kind', [
    'chisquare', 'power_divergence', 'normaltest', 'skewtest', 'kurtosistest',
])
def test_one(kind):
    a = np.random.random(size=30,)
    a_ = da.from_array(a, 3)

    dask_test = getattr(dask.array.stats, kind)
    scipy_test = getattr(scipy.stats, kind)

    result = dask_test(a_)
    expected = scipy_test(a)

    assert isinstance(result, Delayed)
    assert allclose(result.compute(), expected)


@pytest.mark.parametrize('kind, kwargs', [
    ('ttest_ind', {}),
    ('ttest_ind', {'equal_var': False}),
    ('ttest_1samp', {}),
    ('ttest_rel', {}),
    ('chisquare', {}),
    ('power_divergence', {}),
    ('power_divergence', {'lambda_': 0}),
    ('power_divergence', {'lambda_': -1}),
    ('power_divergence', {'lambda_': 'neyman'}),
])
def test_two(kind, kwargs):
    a = np.random.random(size=30,)
    b = np.random.random(size=30,)
    a_ = da.from_array(a, 3)
    b_ = da.from_array(b, 3)

    dask_test = getattr(dask.array.stats, kind)
    scipy_test = getattr(scipy.stats, kind)

    with pytest.warns(None):  # maybe overflow warning (powrer_divergence)
        result = dask_test(a_, b_, **kwargs)
        expected = scipy_test(a, b, **kwargs)

    assert isinstance(result, Delayed)
    assert allclose(result.compute(), expected)
    # fails occasionally. shouldn't this be exact?
    # assert dask.compute(*result) == expected


@pytest.mark.parametrize('k', range(5))
def test_moments(k):
    x = np.random.random(size=(30, 2))
    y = da.from_array(x, 3)

    expected = scipy.stats.moment(x, k)
    result = dask.array.stats.moment(y, k)
    assert_eq(result, expected)


def test_anova():
    np_args = [i * np.random.random(size=(30,)) for i in range(4)]
    da_args = [da.from_array(x, chunks=10) for x in np_args]

    result = dask.array.stats.f_oneway(*da_args)
    expected = scipy.stats.f_oneway(*np_args)

    assert allclose(result.compute(), expected)


@pytest.mark.parametrize('func, nargs', [
    (dask.array.stats.ttest_1samp, 2),
    (dask.array.stats.ttest_rel, 2),
    (dask.array.stats.skewtest, 1),
    (dask.array.stats.kurtosis, 1),
    (dask.array.stats.kurtosistest, 1),
    (dask.array.stats.normaltest, 1),
    (dask.array.stats.moment, 1),
])
@pytest.mark.parametrize('nan_policy', ['omit', 'raise'])
def test_nan_raises(func, nargs, nan_policy):
    with pytest.raises(NotImplementedError):
        func(*(None,) * nargs, nan_policy=nan_policy)


def test_power_divergence_invalid():
    a = np.random.random(size=30,)
    a_ = da.from_array(a, 3)

    with pytest.raises(ValueError):
        dask.array.stats.power_divergence(a_, lambda_='wrong')


def test_skew_raises():
    a = da.ones((7,), chunks=(7,))
    with pytest.raises(ValueError) as rec:
        dask.array.stats.skewtest(a)

    assert "7 samples" in str(rec)
