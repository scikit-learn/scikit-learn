# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest

from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from pandas import (Series, DataFrame,
                    date_range, Timestamp, DatetimeIndex, NaT)

from pandas.compat import lrange, range
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal, assert_almost_equal)

import pandas.util.testing as tm

import pandas._libs.index as _index
from pandas._libs import tslib


"""
Also test support for datetime64[ns] in Series / DataFrame
"""


def test_fancy_getitem():
    dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                        end=datetime(2010, 1, 1))

    s = Series(np.arange(len(dti)), index=dti)

    assert s[48] == 48
    assert s['1/2/2009'] == 48
    assert s['2009-1-2'] == 48
    assert s[datetime(2009, 1, 2)] == 48
    assert s[Timestamp(datetime(2009, 1, 2))] == 48
    pytest.raises(KeyError, s.__getitem__, '2009-1-3')

    assert_series_equal(s['3/6/2009':'2009-06-05'],
                        s[datetime(2009, 3, 6):datetime(2009, 6, 5)])


def test_fancy_setitem():
    dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                        end=datetime(2010, 1, 1))

    s = Series(np.arange(len(dti)), index=dti)
    s[48] = -1
    assert s[48] == -1
    s['1/2/2009'] = -2
    assert s[48] == -2
    s['1/2/2009':'2009-06-05'] = -3
    assert (s[48:54] == -3).all()


def test_dti_snap():
    dti = DatetimeIndex(['1/1/2002', '1/2/2002', '1/3/2002', '1/4/2002',
                         '1/5/2002', '1/6/2002', '1/7/2002'], freq='D')

    res = dti.snap(freq='W-MON')
    exp = date_range('12/31/2001', '1/7/2002', freq='w-mon')
    exp = exp.repeat([3, 4])
    assert (res == exp).all()

    res = dti.snap(freq='B')

    exp = date_range('1/1/2002', '1/7/2002', freq='b')
    exp = exp.repeat([1, 1, 1, 2, 2])
    assert (res == exp).all()


def test_dti_reset_index_round_trip():
    dti = DatetimeIndex(start='1/1/2001', end='6/1/2001', freq='D')
    d1 = DataFrame({'v': np.random.rand(len(dti))}, index=dti)
    d2 = d1.reset_index()
    assert d2.dtypes[0] == np.dtype('M8[ns]')
    d3 = d2.set_index('index')
    assert_frame_equal(d1, d3, check_names=False)

    # #2329
    stamp = datetime(2012, 11, 22)
    df = DataFrame([[stamp, 12.1]], columns=['Date', 'Value'])
    df = df.set_index('Date')

    assert df.index[0] == stamp
    assert df.reset_index()['Date'][0] == stamp


def test_series_set_value():
    # #1561

    dates = [datetime(2001, 1, 1), datetime(2001, 1, 2)]
    index = DatetimeIndex(dates)

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        s = Series().set_value(dates[0], 1.)
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        s2 = s.set_value(dates[1], np.nan)

    exp = Series([1., np.nan], index=index)

    assert_series_equal(s2, exp)

    # s = Series(index[:1], index[:1])
    # s2 = s.set_value(dates[1], index[1])
    # assert s2.values.dtype == 'M8[ns]'


@pytest.mark.slow
def test_slice_locs_indexerror():
    times = [datetime(2000, 1, 1) + timedelta(minutes=i * 10)
             for i in range(100000)]
    s = Series(lrange(100000), times)
    s.loc[datetime(1900, 1, 1):datetime(2100, 1, 1)]


def test_slicing_datetimes():
    # GH 7523

    # unique
    df = DataFrame(np.arange(4., dtype='float64'),
                   index=[datetime(2001, 1, i, 10, 00)
                          for i in [1, 2, 3, 4]])
    result = df.loc[datetime(2001, 1, 1, 10):]
    assert_frame_equal(result, df)
    result = df.loc[:datetime(2001, 1, 4, 10)]
    assert_frame_equal(result, df)
    result = df.loc[datetime(2001, 1, 1, 10):datetime(2001, 1, 4, 10)]
    assert_frame_equal(result, df)

    result = df.loc[datetime(2001, 1, 1, 11):]
    expected = df.iloc[1:]
    assert_frame_equal(result, expected)
    result = df.loc['20010101 11':]
    assert_frame_equal(result, expected)

    # duplicates
    df = pd.DataFrame(np.arange(5., dtype='float64'),
                      index=[datetime(2001, 1, i, 10, 00)
                             for i in [1, 2, 2, 3, 4]])

    result = df.loc[datetime(2001, 1, 1, 10):]
    assert_frame_equal(result, df)
    result = df.loc[:datetime(2001, 1, 4, 10)]
    assert_frame_equal(result, df)
    result = df.loc[datetime(2001, 1, 1, 10):datetime(2001, 1, 4, 10)]
    assert_frame_equal(result, df)

    result = df.loc[datetime(2001, 1, 1, 11):]
    expected = df.iloc[1:]
    assert_frame_equal(result, expected)
    result = df.loc['20010101 11':]
    assert_frame_equal(result, expected)


def test_frame_datetime64_duplicated():
    dates = date_range('2010-07-01', end='2010-08-05')

    tst = DataFrame({'symbol': 'AAA', 'date': dates})
    result = tst.duplicated(['date', 'symbol'])
    assert (-result).all()

    tst = DataFrame({'date': dates})
    result = tst.duplicated()
    assert (-result).all()


def test_getitem_setitem_datetime_tz_pytz():
    from pytz import timezone as tz
    from pandas import date_range

    N = 50
    # testing with timezone, GH #2785
    rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
    ts = Series(np.random.randn(N), index=rng)

    # also test Timestamp tz handling, GH #2789
    result = ts.copy()
    result["1990-01-01 09:00:00+00:00"] = 0
    result["1990-01-01 09:00:00+00:00"] = ts[4]
    assert_series_equal(result, ts)

    result = ts.copy()
    result["1990-01-01 03:00:00-06:00"] = 0
    result["1990-01-01 03:00:00-06:00"] = ts[4]
    assert_series_equal(result, ts)

    # repeat with datetimes
    result = ts.copy()
    result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
    result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
    assert_series_equal(result, ts)

    result = ts.copy()

    # comparison dates with datetime MUST be localized!
    date = tz('US/Central').localize(datetime(1990, 1, 1, 3))
    result[date] = 0
    result[date] = ts[4]
    assert_series_equal(result, ts)


def test_getitem_setitem_datetime_tz_dateutil():
    from dateutil.tz import tzutc
    from pandas._libs.tslibs.timezones import dateutil_gettz as gettz

    tz = lambda x: tzutc() if x == 'UTC' else gettz(
        x)  # handle special case for utc in dateutil

    from pandas import date_range

    N = 50

    # testing with timezone, GH #2785
    rng = date_range('1/1/1990', periods=N, freq='H',
                     tz='America/New_York')
    ts = Series(np.random.randn(N), index=rng)

    # also test Timestamp tz handling, GH #2789
    result = ts.copy()
    result["1990-01-01 09:00:00+00:00"] = 0
    result["1990-01-01 09:00:00+00:00"] = ts[4]
    assert_series_equal(result, ts)

    result = ts.copy()
    result["1990-01-01 03:00:00-06:00"] = 0
    result["1990-01-01 03:00:00-06:00"] = ts[4]
    assert_series_equal(result, ts)

    # repeat with datetimes
    result = ts.copy()
    result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
    result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
    assert_series_equal(result, ts)

    result = ts.copy()
    result[datetime(1990, 1, 1, 3, tzinfo=tz('America/Chicago'))] = 0
    result[datetime(1990, 1, 1, 3, tzinfo=tz('America/Chicago'))] = ts[4]
    assert_series_equal(result, ts)


def test_getitem_setitem_datetimeindex():
    N = 50
    # testing with timezone, GH #2785
    rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
    ts = Series(np.random.randn(N), index=rng)

    result = ts["1990-01-01 04:00:00"]
    expected = ts[4]
    assert result == expected

    result = ts.copy()
    result["1990-01-01 04:00:00"] = 0
    result["1990-01-01 04:00:00"] = ts[4]
    assert_series_equal(result, ts)

    result = ts["1990-01-01 04:00:00":"1990-01-01 07:00:00"]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts.copy()
    result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = 0
    result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = ts[4:8]
    assert_series_equal(result, ts)

    lb = "1990-01-01 04:00:00"
    rb = "1990-01-01 07:00:00"
    # GH#18435 strings get a pass from tzawareness compat
    result = ts[(ts.index >= lb) & (ts.index <= rb)]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    lb = "1990-01-01 04:00:00-0500"
    rb = "1990-01-01 07:00:00-0500"
    result = ts[(ts.index >= lb) & (ts.index <= rb)]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    # repeat all the above with naive datetimes
    result = ts[datetime(1990, 1, 1, 4)]
    expected = ts[4]
    assert result == expected

    result = ts.copy()
    result[datetime(1990, 1, 1, 4)] = 0
    result[datetime(1990, 1, 1, 4)] = ts[4]
    assert_series_equal(result, ts)

    result = ts[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts.copy()
    result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = 0
    result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = ts[4:8]
    assert_series_equal(result, ts)

    lb = datetime(1990, 1, 1, 4)
    rb = datetime(1990, 1, 1, 7)
    with pytest.raises(TypeError):
        # tznaive vs tzaware comparison is invalid
        # see GH#18376, GH#18162
        ts[(ts.index >= lb) & (ts.index <= rb)]

    lb = pd.Timestamp(datetime(1990, 1, 1, 4)).tz_localize(rng.tzinfo)
    rb = pd.Timestamp(datetime(1990, 1, 1, 7)).tz_localize(rng.tzinfo)
    result = ts[(ts.index >= lb) & (ts.index <= rb)]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts[ts.index[4]]
    expected = ts[4]
    assert result == expected

    result = ts[ts.index[4:8]]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts.copy()
    result[ts.index[4:8]] = 0
    result[4:8] = ts[4:8]
    assert_series_equal(result, ts)

    # also test partial date slicing
    result = ts["1990-01-02"]
    expected = ts[24:48]
    assert_series_equal(result, expected)

    result = ts.copy()
    result["1990-01-02"] = 0
    result["1990-01-02"] = ts[24:48]
    assert_series_equal(result, ts)


def test_getitem_setitem_periodindex():
    from pandas import period_range

    N = 50
    rng = period_range('1/1/1990', periods=N, freq='H')
    ts = Series(np.random.randn(N), index=rng)

    result = ts["1990-01-01 04"]
    expected = ts[4]
    assert result == expected

    result = ts.copy()
    result["1990-01-01 04"] = 0
    result["1990-01-01 04"] = ts[4]
    assert_series_equal(result, ts)

    result = ts["1990-01-01 04":"1990-01-01 07"]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts.copy()
    result["1990-01-01 04":"1990-01-01 07"] = 0
    result["1990-01-01 04":"1990-01-01 07"] = ts[4:8]
    assert_series_equal(result, ts)

    lb = "1990-01-01 04"
    rb = "1990-01-01 07"
    result = ts[(ts.index >= lb) & (ts.index <= rb)]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    # GH 2782
    result = ts[ts.index[4]]
    expected = ts[4]
    assert result == expected

    result = ts[ts.index[4:8]]
    expected = ts[4:8]
    assert_series_equal(result, expected)

    result = ts.copy()
    result[ts.index[4:8]] = 0
    result[4:8] = ts[4:8]
    assert_series_equal(result, ts)


def test_getitem_median_slice_bug():
    index = date_range('20090415', '20090519', freq='2B')
    s = Series(np.random.randn(13), index=index)

    indexer = [slice(6, 7, None)]
    result = s[indexer]
    expected = s[indexer[0]]
    assert_series_equal(result, expected)


def test_datetime_indexing():
    from pandas import date_range

    index = date_range('1/1/2000', '1/7/2000')
    index = index.repeat(3)

    s = Series(len(index), index=index)
    stamp = Timestamp('1/8/2000')

    pytest.raises(KeyError, s.__getitem__, stamp)
    s[stamp] = 0
    assert s[stamp] == 0

    # not monotonic
    s = Series(len(index), index=index)
    s = s[::-1]

    pytest.raises(KeyError, s.__getitem__, stamp)
    s[stamp] = 0
    assert s[stamp] == 0


"""
test duplicates in time series
"""


@pytest.fixture(scope='module')
def dups():
    dates = [datetime(2000, 1, 2), datetime(2000, 1, 2),
             datetime(2000, 1, 2), datetime(2000, 1, 3),
             datetime(2000, 1, 3), datetime(2000, 1, 3),
             datetime(2000, 1, 4), datetime(2000, 1, 4),
             datetime(2000, 1, 4), datetime(2000, 1, 5)]

    return Series(np.random.randn(len(dates)), index=dates)


def test_constructor(dups):
    assert isinstance(dups, Series)
    assert isinstance(dups.index, DatetimeIndex)


def test_is_unique_monotonic(dups):
    assert not dups.index.is_unique


def test_index_unique(dups):
    uniques = dups.index.unique()
    expected = DatetimeIndex([datetime(2000, 1, 2), datetime(2000, 1, 3),
                              datetime(2000, 1, 4), datetime(2000, 1, 5)])
    assert uniques.dtype == 'M8[ns]'  # sanity
    tm.assert_index_equal(uniques, expected)
    assert dups.index.nunique() == 4

    # #2563
    assert isinstance(uniques, DatetimeIndex)

    dups_local = dups.index.tz_localize('US/Eastern')
    dups_local.name = 'foo'
    result = dups_local.unique()
    expected = DatetimeIndex(expected, name='foo')
    expected = expected.tz_localize('US/Eastern')
    assert result.tz is not None
    assert result.name == 'foo'
    tm.assert_index_equal(result, expected)

    # NaT, note this is excluded
    arr = [1370745748 + t for t in range(20)] + [tslib.iNaT]
    idx = DatetimeIndex(arr * 3)
    tm.assert_index_equal(idx.unique(), DatetimeIndex(arr))
    assert idx.nunique() == 20
    assert idx.nunique(dropna=False) == 21

    arr = [Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t)
           for t in range(20)] + [NaT]
    idx = DatetimeIndex(arr * 3)
    tm.assert_index_equal(idx.unique(), DatetimeIndex(arr))
    assert idx.nunique() == 20
    assert idx.nunique(dropna=False) == 21


def test_index_dupes_contains():
    d = datetime(2011, 12, 5, 20, 30)
    ix = DatetimeIndex([d, d])
    assert d in ix


def test_duplicate_dates_indexing(dups):
    ts = dups

    uniques = ts.index.unique()
    for date in uniques:
        result = ts[date]

        mask = ts.index == date
        total = (ts.index == date).sum()
        expected = ts[mask]
        if total > 1:
            assert_series_equal(result, expected)
        else:
            assert_almost_equal(result, expected[0])

        cp = ts.copy()
        cp[date] = 0
        expected = Series(np.where(mask, 0, ts), index=ts.index)
        assert_series_equal(cp, expected)

    pytest.raises(KeyError, ts.__getitem__, datetime(2000, 1, 6))

    # new index
    ts[datetime(2000, 1, 6)] = 0
    assert ts[datetime(2000, 1, 6)] == 0


def test_range_slice():
    idx = DatetimeIndex(['1/1/2000', '1/2/2000', '1/2/2000', '1/3/2000',
                         '1/4/2000'])

    ts = Series(np.random.randn(len(idx)), index=idx)

    result = ts['1/2/2000':]
    expected = ts[1:]
    assert_series_equal(result, expected)

    result = ts['1/2/2000':'1/3/2000']
    expected = ts[1:4]
    assert_series_equal(result, expected)


def test_groupby_average_dup_values(dups):
    result = dups.groupby(level=0).mean()
    expected = dups.groupby(dups.index).mean()
    assert_series_equal(result, expected)


def test_indexing_over_size_cutoff():
    import datetime
    # #1821

    old_cutoff = _index._SIZE_CUTOFF
    try:
        _index._SIZE_CUTOFF = 1000

        # create large list of non periodic datetime
        dates = []
        sec = datetime.timedelta(seconds=1)
        half_sec = datetime.timedelta(microseconds=500000)
        d = datetime.datetime(2011, 12, 5, 20, 30)
        n = 1100
        for i in range(n):
            dates.append(d)
            dates.append(d + sec)
            dates.append(d + sec + half_sec)
            dates.append(d + sec + sec + half_sec)
            d += 3 * sec

        # duplicate some values in the list
        duplicate_positions = np.random.randint(0, len(dates) - 1, 20)
        for p in duplicate_positions:
            dates[p + 1] = dates[p]

        df = DataFrame(np.random.randn(len(dates), 4),
                       index=dates,
                       columns=list('ABCD'))

        pos = n * 3
        timestamp = df.index[pos]
        assert timestamp in df.index

        # it works!
        df.loc[timestamp]
        assert len(df.loc[[timestamp]]) > 0
    finally:
        _index._SIZE_CUTOFF = old_cutoff


def test_indexing_unordered():
    # GH 2437
    rng = date_range(start='2011-01-01', end='2011-01-15')
    ts = Series(np.random.rand(len(rng)), index=rng)
    ts2 = pd.concat([ts[0:4], ts[-4:], ts[4:-4]])

    for t in ts.index:
        # TODO: unused?
        s = str(t)  # noqa

        expected = ts[t]
        result = ts2[t]
        assert expected == result

    # GH 3448 (ranges)
    def compare(slobj):
        result = ts2[slobj].copy()
        result = result.sort_index()
        expected = ts[slobj]
        assert_series_equal(result, expected)

    compare(slice('2011-01-01', '2011-01-15'))
    compare(slice('2010-12-30', '2011-01-15'))
    compare(slice('2011-01-01', '2011-01-16'))

    # partial ranges
    compare(slice('2011-01-01', '2011-01-6'))
    compare(slice('2011-01-06', '2011-01-8'))
    compare(slice('2011-01-06', '2011-01-12'))

    # single values
    result = ts2['2011'].sort_index()
    expected = ts['2011']
    assert_series_equal(result, expected)

    # diff freq
    rng = date_range(datetime(2005, 1, 1), periods=20, freq='M')
    ts = Series(np.arange(len(rng)), index=rng)
    ts = ts.take(np.random.permutation(20))

    result = ts['2005']
    for t in result.index:
        assert t.year == 2005


def test_indexing():
    idx = date_range("2001-1-1", periods=20, freq='M')
    ts = Series(np.random.rand(len(idx)), index=idx)

    # getting

    # GH 3070, make sure semantics work on Series/Frame
    expected = ts['2001']
    expected.name = 'A'

    df = DataFrame(dict(A=ts))
    result = df['2001']['A']
    assert_series_equal(expected, result)

    # setting
    ts['2001'] = 1
    expected = ts['2001']
    expected.name = 'A'

    df.loc['2001', 'A'] = 1

    result = df['2001']['A']
    assert_series_equal(expected, result)

    # GH3546 (not including times on the last day)
    idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:00',
                     freq='H')
    ts = Series(lrange(len(idx)), index=idx)
    expected = ts['2013-05']
    assert_series_equal(expected, ts)

    idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:59',
                     freq='S')
    ts = Series(lrange(len(idx)), index=idx)
    expected = ts['2013-05']
    assert_series_equal(expected, ts)

    idx = [Timestamp('2013-05-31 00:00'),
           Timestamp(datetime(2013, 5, 31, 23, 59, 59, 999999))]
    ts = Series(lrange(len(idx)), index=idx)
    expected = ts['2013']
    assert_series_equal(expected, ts)

    # GH14826, indexing with a seconds resolution string / datetime object
    df = DataFrame(np.random.rand(5, 5),
                   columns=['open', 'high', 'low', 'close', 'volume'],
                   index=date_range('2012-01-02 18:01:00',
                                    periods=5, tz='US/Central', freq='s'))
    expected = df.loc[[df.index[2]]]

    # this is a single date, so will raise
    pytest.raises(KeyError, df.__getitem__, '2012-01-02 18:01:02', )
    pytest.raises(KeyError, df.__getitem__, df.index[2], )


"""
test NaT support
"""


def test_set_none_nan():
    series = Series(date_range('1/1/2000', periods=10))
    series[3] = None
    assert series[3] is NaT

    series[3:5] = None
    assert series[4] is NaT

    series[5] = np.nan
    assert series[5] is NaT

    series[5:7] = np.nan
    assert series[6] is NaT


def test_nat_operations():
    # GH 8617
    s = Series([0, pd.NaT], dtype='m8[ns]')
    exp = s[0]
    assert s.median() == exp
    assert s.min() == exp
    assert s.max() == exp


@pytest.mark.parametrize('method', ["round", "floor", "ceil"])
@pytest.mark.parametrize('freq', ["s", "5s", "min", "5min", "h", "5h"])
def test_round_nat(method, freq):
    # GH14940
    s = Series([pd.NaT])
    expected = Series(pd.NaT)
    round_method = getattr(s.dt, method)
    assert_series_equal(round_method(freq), expected)
