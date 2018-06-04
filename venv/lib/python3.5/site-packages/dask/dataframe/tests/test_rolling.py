import pandas as pd
import pytest
import numpy as np

import dask.dataframe as dd
from dask.dataframe.utils import assert_eq, PANDAS_VERSION

N = 40
df = pd.DataFrame({'a': np.random.randn(N).cumsum(),
                   'b': np.random.randint(100, size=(N,)),
                   'c': np.random.randint(100, size=(N,)),
                   'd': np.random.randint(100, size=(N,)),
                   'e': np.random.randint(100, size=(N,))})
ddf = dd.from_pandas(df, 3)

idx = (pd.date_range('2016-01-01', freq='3s', periods=100) |
       pd.date_range('2016-01-01', freq='5s', periods=100))[:N]

ts = pd.DataFrame({'a': np.random.randn(N).cumsum(),
                   'b': np.random.randint(100, size=(N,)),
                   'c': np.random.randint(100, size=(N,)),
                   'd': np.random.randint(100, size=(N,)),
                   'e': np.random.randint(100, size=(N,))},
                  index=idx)
dts = dd.from_pandas(ts, 3)


def shifted_sum(df, before, after, c=0):
    a = df.shift(before)
    b = df.shift(-after)
    return df + a + b + c


def ts_shifted_sum(df, before, after, c=0):
    a = df.shift(before.seconds)
    b = df.shift(-after.seconds)
    return df + a + b + c


@pytest.mark.parametrize('npartitions', [1, 4])
def test_map_overlap(npartitions):
    ddf = dd.from_pandas(df, npartitions)
    for before, after in [(0, 3), (3, 0), (3, 3), (0, 0)]:
        # DataFrame
        res = ddf.map_overlap(shifted_sum, before, after, before, after, c=2)
        sol = shifted_sum(df, before, after, c=2)
        assert_eq(res, sol)

        # Series
        res = ddf.b.map_overlap(shifted_sum, before, after, before, after, c=2)
        sol = shifted_sum(df.b, before, after, c=2)
        assert_eq(res, sol)


def test_map_partitions_names():
    npartitions = 3
    ddf = dd.from_pandas(df, npartitions)

    res = ddf.map_overlap(shifted_sum, 0, 3, 0, 3, c=2)
    res2 = ddf.map_overlap(shifted_sum, 0, 3, 0, 3, c=2)
    assert set(res.dask) == set(res2.dask)

    res3 = ddf.map_overlap(shifted_sum, 0, 3, 0, 3, c=3)
    assert res3._name != res._name
    # Difference is just the final map
    diff = set(res3.dask).difference(res.dask)
    assert len(diff) == npartitions

    res4 = ddf.map_overlap(shifted_sum, 3, 0, 0, 3, c=2)
    assert res4._name != res._name


def test_map_partitions_errors():
    # Non-integer
    with pytest.raises(ValueError):
        ddf.map_overlap(shifted_sum, 0.5, 3, 0, 2, c=2)

    # Negative
    with pytest.raises(ValueError):
        ddf.map_overlap(shifted_sum, 0, -5, 0, 2, c=2)

    # Partition size < window size
    with pytest.raises(NotImplementedError):
        ddf.map_overlap(shifted_sum, 0, 100, 0, 100, c=2).compute()

    # Offset with non-datetime
    with pytest.raises(TypeError):
        ddf.map_overlap(shifted_sum, pd.Timedelta('1s'), pd.Timedelta('1s'),
                        0, 2, c=2)


def mad(x):
    return np.fabs(x - x.mean()).mean()


rolling_method_args_check_less_precise = [
    ('count', (), False),
    ('sum', (), False),
    ('mean', (), False),
    ('median', (), False),
    ('min', (), False),
    ('max', (), False),
    ('std', (), False),
    ('var', (), False),
    ('skew', (), True),   # here and elsewhere, results for kurt and skew are
    ('kurt', (), True),   # checked with check_less_precise=True so that we are
                          # only looking at 3ish decimal places for the equality check
                          # rather than 5ish. I have encountered a case where a test
                          # seems to have failed due to numerical problems with kurt.
                          # So far, I am only weakening the check for kurt and skew,
                          # as they involve third degree powers and higher
    ('quantile', (.38,), False),
    ('apply', (mad,), False),
]


@pytest.mark.parametrize('method,args,check_less_precise',
                         rolling_method_args_check_less_precise)
@pytest.mark.parametrize('window', [1, 2, 4, 5])
@pytest.mark.parametrize('center', [True, False])
def test_rolling_methods(method, args, window, center, check_less_precise):
    # DataFrame
    prolling = df.rolling(window, center=center)
    drolling = ddf.rolling(window, center=center)
    if method == 'apply' and PANDAS_VERSION >= '0.23.0':
        kwargs = {'raw': False}
    else:
        kwargs = {}
    assert_eq(getattr(prolling, method)(*args, **kwargs),
              getattr(drolling, method)(*args, **kwargs),
              check_less_precise=check_less_precise)

    # Series
    prolling = df.a.rolling(window, center=center)
    drolling = ddf.a.rolling(window, center=center)
    assert_eq(getattr(prolling, method)(*args, **kwargs),
              getattr(drolling, method)(*args, **kwargs),
              check_less_precise=check_less_precise)


@pytest.mark.skipif(PANDAS_VERSION >= '0.23.0', reason="Raw is allowed.")
def test_rolling_raw_pandas_lt_0230_raises():
    with pytest.raises(TypeError):
        df.rolling(2).apply(mad, raw=True)


def test_rolling_raises():
    df = pd.DataFrame({'a': np.random.randn(25).cumsum(),
                       'b': np.random.randint(100, size=(25,))})
    ddf = dd.from_pandas(df, 3)
    pytest.raises(ValueError, lambda: ddf.rolling(1.5))
    pytest.raises(ValueError, lambda: ddf.rolling(-1))
    pytest.raises(ValueError, lambda: ddf.rolling(3, min_periods=1.2))
    pytest.raises(ValueError, lambda: ddf.rolling(3, min_periods=-2))
    pytest.raises(ValueError, lambda: ddf.rolling(3, axis=10))
    pytest.raises(ValueError, lambda: ddf.rolling(3, axis='coulombs'))
    pytest.raises(NotImplementedError, lambda: ddf.rolling(100).mean().compute())


def test_rolling_names():
    df = pd.DataFrame({'a': [1, 2, 3],
                       'b': [4, 5, 6]})
    a = dd.from_pandas(df, npartitions=2)
    assert sorted(a.rolling(2).sum().dask) == sorted(a.rolling(2).sum().dask)


def test_rolling_axis():
    df = pd.DataFrame(np.random.randn(20, 16))
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(df.rolling(3, axis=0).mean(), ddf.rolling(3, axis=0).mean())
    assert_eq(df.rolling(3, axis=1).mean(), ddf.rolling(3, axis=1).mean())
    assert_eq(df.rolling(3, min_periods=1, axis=1).mean(),
              ddf.rolling(3, min_periods=1, axis=1).mean())
    assert_eq(df.rolling(3, axis='columns').mean(),
              ddf.rolling(3, axis='columns').mean())
    assert_eq(df.rolling(3, axis='rows').mean(),
              ddf.rolling(3, axis='rows').mean())

    s = df[3]
    ds = ddf[3]
    assert_eq(s.rolling(5, axis=0).std(), ds.rolling(5, axis=0).std())


def test_rolling_partition_size():
    df = pd.DataFrame(np.random.randn(50, 2))
    ddf = dd.from_pandas(df, npartitions=5)

    for obj, dobj in [(df, ddf), (df[0], ddf[0])]:
        assert_eq(obj.rolling(10).mean(), dobj.rolling(10).mean())
        assert_eq(obj.rolling(11).mean(), dobj.rolling(11).mean())
        with pytest.raises(NotImplementedError):
            dobj.rolling(12).mean().compute()


def test_rolling_repr():
    ddf = dd.from_pandas(pd.DataFrame([10] * 30), npartitions=3)
    assert repr(ddf.rolling(4)) == 'Rolling [window=4,center=False,axis=0]'


def test_time_rolling_repr():
    assert repr(dts.rolling('4s')) == (
        'Rolling [window=4000000000,center=False,win_type=freq,axis=0]')


def test_time_rolling_constructor():
    result = dts.rolling('4s')
    assert result.window == '4s'
    assert result.min_periods is None
    assert result.win_type is None

    assert result._win_type == 'freq'
    assert result._window == 4000000000  # ns
    assert result._min_periods == 1


@pytest.mark.parametrize('method,args,check_less_precise',
                         rolling_method_args_check_less_precise)
@pytest.mark.parametrize('window', ['1S', '2S', '3S', pd.offsets.Second(5)])
def test_time_rolling_methods(method, args, window, check_less_precise):
    # DataFrame
    if method == 'apply' and PANDAS_VERSION >= '0.23.0':
        kwargs = {"raw": False}
    else:
        kwargs = {}
    prolling = ts.rolling(window)
    drolling = dts.rolling(window)
    assert_eq(getattr(prolling, method)(*args, **kwargs),
              getattr(drolling, method)(*args, **kwargs),
              check_less_precise=check_less_precise)

    # Series
    prolling = ts.a.rolling(window)
    drolling = dts.a.rolling(window)
    assert_eq(getattr(prolling, method)(*args, **kwargs),
              getattr(drolling, method)(*args, **kwargs),
              check_less_precise=check_less_precise)


@pytest.mark.parametrize('window', [pd.Timedelta('31s'), pd.Timedelta('1M')])
def test_time_rolling_window_too_large(window):
    with pytest.raises(ValueError):
        dts.map_overlap(ts_shifted_sum, window, window, window, window, c=2)


@pytest.mark.parametrize('before, after', [
    ('6s', '6s'),
    ('2s', '2s'),
    ('6s', '2s'),
])
def test_time_rolling(before, after):
    window = before
    before = pd.Timedelta(before)
    after = pd.Timedelta(after)
    result = dts.map_overlap(lambda x: x.rolling(window).count(), before, after)
    expected = dts.compute().rolling(window).count()
    assert_eq(result, expected)
