import pandas.util.testing as tm
import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe.utils import assert_eq


def test_make_timeseries():
    df = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                 freq='2D', partition_freq='6M')

    assert df.divisions[0] == pd.Timestamp('2000-01-31', freq='6M')
    assert df.divisions[-1] == pd.Timestamp('2014-07-31', freq='6M')
    tm.assert_index_equal(df.columns, pd.Index(['A', 'B', 'C']))
    assert df['A'].head().dtype == float
    assert df['B'].head().dtype == int
    assert df['C'].head().dtype == object
    assert df.divisions == tuple(pd.DatetimeIndex(start='2000', end='2015',
                                                  freq='6M'))

    tm.assert_frame_equal(df.head(), df.head())

    a = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                freq='2D', partition_freq='6M', seed=123)
    b = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                freq='2D', partition_freq='6M', seed=123)
    c = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                freq='2D', partition_freq='6M', seed=456)
    d = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                freq='2D', partition_freq='3M', seed=123)
    e = dd.demo.make_timeseries('2000', '2015', {'A': float, 'B': int, 'C': str},
                                freq='1D', partition_freq='6M', seed=123)
    tm.assert_frame_equal(a.head(), b.head())
    assert not (a.head(10) == c.head(10)).all().all()
    assert a._name == b._name
    assert a._name != c._name
    assert a._name != d._name
    assert a._name != e._name


def test_make_timeseries_no_args():
    df = dd.demo.make_timeseries()
    assert 1 < df.npartitions < 1000
    assert len(df.columns) > 1
    assert len(set(df.dtypes)) > 1


def test_no_overlaps():
    df = dd.demo.make_timeseries('2000', '2001', {'A': float},
                                 freq='3H', partition_freq='3M')

    assert all(df.get_partition(i).index.max().compute() <
               df.get_partition(i + 1).index.min().compute()
               for i in range(df.npartitions - 2))


@pytest.mark.xfail
@pytest.mark.network
def test_daily_stock():
    pytest.importorskip('pandas_datareader')
    df = dd.demo.daily_stock('GOOG', start='2010-01-01', stop='2010-01-30', freq='1h')
    assert isinstance(df, dd.DataFrame)
    assert 10 < df.npartitions < 31
    assert_eq(df, df)
