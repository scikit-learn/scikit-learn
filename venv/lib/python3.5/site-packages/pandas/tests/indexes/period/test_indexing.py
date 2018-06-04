from datetime import datetime, timedelta

import pytest

import numpy as np
import pandas as pd
from pandas.util import testing as tm
from pandas.compat import lrange
from pandas._libs import tslibs
from pandas import (PeriodIndex, Series, DatetimeIndex,
                    period_range, Period, notna)
from pandas._libs.tslibs import period as libperiod


class TestGetItem(object):

    def test_getitem(self):
        idx1 = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')

        for idx in [idx1]:
            result = idx[0]
            assert result == pd.Period('2011-01-01', freq='D')

            result = idx[-1]
            assert result == pd.Period('2011-01-31', freq='D')

            result = idx[0:5]
            expected = pd.period_range('2011-01-01', '2011-01-05', freq='D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx[0:10:2]
            expected = pd.PeriodIndex(['2011-01-01', '2011-01-03',
                                       '2011-01-05',
                                       '2011-01-07', '2011-01-09'],
                                      freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx[-20:-5:3]
            expected = pd.PeriodIndex(['2011-01-12', '2011-01-15',
                                       '2011-01-18',
                                       '2011-01-21', '2011-01-24'],
                                      freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx[4::-1]
            expected = PeriodIndex(['2011-01-05', '2011-01-04', '2011-01-03',
                                    '2011-01-02', '2011-01-01'],
                                   freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

    def test_getitem_index(self):
        idx = period_range('2007-01', periods=10, freq='M', name='x')

        result = idx[[1, 3, 5]]
        exp = pd.PeriodIndex(['2007-02', '2007-04', '2007-06'],
                             freq='M', name='x')
        tm.assert_index_equal(result, exp)

        result = idx[[True, True, False, False, False,
                      True, True, False, False, False]]
        exp = pd.PeriodIndex(['2007-01', '2007-02', '2007-06', '2007-07'],
                             freq='M', name='x')
        tm.assert_index_equal(result, exp)

    def test_getitem_partial(self):
        rng = period_range('2007-01', periods=50, freq='M')
        ts = Series(np.random.randn(len(rng)), rng)

        pytest.raises(KeyError, ts.__getitem__, '2006')

        result = ts['2008']
        assert (result.index.year == 2008).all()

        result = ts['2008':'2009']
        assert len(result) == 24

        result = ts['2008-1':'2009-12']
        assert len(result) == 24

        result = ts['2008Q1':'2009Q4']
        assert len(result) == 24

        result = ts[:'2009']
        assert len(result) == 36

        result = ts['2009':]
        assert len(result) == 50 - 24

        exp = result
        result = ts[24:]
        tm.assert_series_equal(exp, result)

        ts = ts[10:].append(ts[10:])
        tm.assert_raises_regex(KeyError,
                               "left slice bound for non-unique "
                               "label: '2008'",
                               ts.__getitem__, slice('2008', '2009'))

    def test_getitem_datetime(self):
        rng = period_range(start='2012-01-01', periods=10, freq='W-MON')
        ts = Series(lrange(len(rng)), index=rng)

        dt1 = datetime(2011, 10, 2)
        dt4 = datetime(2012, 4, 20)

        rs = ts[dt1:dt4]
        tm.assert_series_equal(rs, ts)

    def test_getitem_nat(self):
        idx = pd.PeriodIndex(['2011-01', 'NaT', '2011-02'], freq='M')
        assert idx[0] == pd.Period('2011-01', freq='M')
        assert idx[1] is pd.NaT

        s = pd.Series([0, 1, 2], index=idx)
        assert s[pd.NaT] == 1

        s = pd.Series(idx, index=idx)
        assert (s[pd.Period('2011-01', freq='M')] ==
                pd.Period('2011-01', freq='M'))
        assert s[pd.NaT] is pd.NaT

    def test_getitem_list_periods(self):
        # GH 7710
        rng = period_range(start='2012-01-01', periods=10, freq='D')
        ts = Series(lrange(len(rng)), index=rng)
        exp = ts.iloc[[1]]
        tm.assert_series_equal(ts[[Period('2012-01-02', freq='D')]], exp)

    def test_getitem_seconds(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01 09:00:00', freq='S',
                             periods=4000)
        pidx = PeriodIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:
                # GH7116
                # these show deprecations as we are trying
                # to slice with non-integer indexers
                # with pytest.raises(IndexError):
                #    idx[v]
                continue

            s = Series(np.random.rand(len(idx)), index=idx)
            tm.assert_series_equal(s['2013/01/01 10:00'], s[3600:3660])
            tm.assert_series_equal(s['2013/01/01 9H'], s[:3600])
            for d in ['2013/01/01', '2013/01', '2013']:
                tm.assert_series_equal(s[d], s)

    def test_getitem_day(self):
        # GH 6716
        # Confirm DatetimeIndex and PeriodIndex works identically
        didx = DatetimeIndex(start='2013/01/01', freq='D', periods=400)
        pidx = PeriodIndex(start='2013/01/01', freq='D', periods=400)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:

                # GH7116
                # these show deprecations as we are trying
                # to slice with non-integer indexers
                # with pytest.raises(IndexError):
                #    idx[v]
                continue

            s = Series(np.random.rand(len(idx)), index=idx)
            tm.assert_series_equal(s['2013/01'], s[0:31])
            tm.assert_series_equal(s['2013/02'], s[31:59])
            tm.assert_series_equal(s['2014'], s[365:])

            invalid = ['2013/02/01 9H', '2013/02/01 09:00']
            for v in invalid:
                with pytest.raises(KeyError):
                    s[v]


class TestWhere(object):
    @pytest.mark.parametrize('klass', [list, tuple, np.array, Series])
    def test_where(self, klass):
        i = period_range('20130101', periods=5, freq='D')
        cond = [True] * len(i)
        expected = i
        result = i.where(klass(cond))
        tm.assert_index_equal(result, expected)

        cond = [False] + [True] * (len(i) - 1)
        expected = PeriodIndex([pd.NaT] + i[1:].tolist(), freq='D')
        result = i.where(klass(cond))
        tm.assert_index_equal(result, expected)

    def test_where_other(self):
        i = period_range('20130101', periods=5, freq='D')
        for arr in [np.nan, pd.NaT]:
            result = i.where(notna(i), other=np.nan)
            expected = i
            tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = pd.PeriodIndex([pd.NaT, pd.NaT] + i[2:].tolist(),
                            freq='D')
        result = i.where(notna(i2), i2)
        tm.assert_index_equal(result, i2)

        i2 = i.copy()
        i2 = pd.PeriodIndex([pd.NaT, pd.NaT] + i[2:].tolist(),
                            freq='D')
        result = i.where(notna(i2), i2.values)
        tm.assert_index_equal(result, i2)


class TestTake(object):
    def test_take(self):
        # GH#10295
        idx1 = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')

        for idx in [idx1]:
            result = idx.take([0])
            assert result == pd.Period('2011-01-01', freq='D')

            result = idx.take([5])
            assert result == pd.Period('2011-01-06', freq='D')

            result = idx.take([0, 1, 2])
            expected = pd.period_range('2011-01-01', '2011-01-03', freq='D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == 'D'
            assert result.freq == expected.freq

            result = idx.take([0, 2, 4])
            expected = pd.PeriodIndex(['2011-01-01', '2011-01-03',
                                       '2011-01-05'], freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx.take([7, 4, 1])
            expected = pd.PeriodIndex(['2011-01-08', '2011-01-05',
                                       '2011-01-02'],
                                      freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx.take([3, 2, 5])
            expected = PeriodIndex(['2011-01-04', '2011-01-03', '2011-01-06'],
                                   freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

            result = idx.take([-3, 2, 5])
            expected = PeriodIndex(['2011-01-29', '2011-01-03', '2011-01-06'],
                                   freq='D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq
            assert result.freq == 'D'

    def test_take_misc(self):
        index = PeriodIndex(start='1/1/10', end='12/31/12', freq='D',
                            name='idx')
        expected = PeriodIndex([datetime(2010, 1, 6), datetime(2010, 1, 7),
                                datetime(2010, 1, 9), datetime(2010, 1, 13)],
                               freq='D', name='idx')

        taken1 = index.take([5, 6, 8, 12])
        taken2 = index[[5, 6, 8, 12]]

        for taken in [taken1, taken2]:
            tm.assert_index_equal(taken, expected)
            assert isinstance(taken, PeriodIndex)
            assert taken.freq == index.freq
            assert taken.name == expected.name

    def test_take_fill_value(self):
        # GH#12631
        idx = pd.PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                             name='xxx', freq='D')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with pytest.raises(IndexError):
            idx.take(np.array([1, -5]))


class TestIndexing(object):

    def test_get_loc_msg(self):
        idx = period_range('2000-1-1', freq='A', periods=10)
        bad_period = Period('2012', 'A')
        pytest.raises(KeyError, idx.get_loc, bad_period)

        try:
            idx.get_loc(bad_period)
        except KeyError as inst:
            assert inst.args[0] == bad_period

    def test_get_loc_nat(self):
        didx = DatetimeIndex(['2011-01-01', 'NaT', '2011-01-03'])
        pidx = PeriodIndex(['2011-01-01', 'NaT', '2011-01-03'], freq='M')

        # check DatetimeIndex compat
        for idx in [didx, pidx]:
            assert idx.get_loc(pd.NaT) == 1
            assert idx.get_loc(None) == 1
            assert idx.get_loc(float('nan')) == 1
            assert idx.get_loc(np.nan) == 1

    def test_get_loc(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')

        # get the location of p1/p2 from
        # monotonic increasing PeriodIndex with non-duplicate
        idx0 = pd.PeriodIndex([p0, p1, p2])
        expected_idx1_p1 = 1
        expected_idx1_p2 = 2

        assert idx0.get_loc(p1) == expected_idx1_p1
        assert idx0.get_loc(str(p1)) == expected_idx1_p1
        assert idx0.get_loc(p2) == expected_idx1_p2
        assert idx0.get_loc(str(p2)) == expected_idx1_p2

        pytest.raises(tslibs.parsing.DateParseError, idx0.get_loc, 'foo')
        pytest.raises(KeyError, idx0.get_loc, 1.1)
        pytest.raises(TypeError, idx0.get_loc, idx0)

        # get the location of p1/p2 from
        # monotonic increasing PeriodIndex with duplicate
        idx1 = pd.PeriodIndex([p1, p1, p2])
        expected_idx1_p1 = slice(0, 2)
        expected_idx1_p2 = 2

        assert idx1.get_loc(p1) == expected_idx1_p1
        assert idx1.get_loc(str(p1)) == expected_idx1_p1
        assert idx1.get_loc(p2) == expected_idx1_p2
        assert idx1.get_loc(str(p2)) == expected_idx1_p2

        pytest.raises(tslibs.parsing.DateParseError, idx1.get_loc, 'foo')
        pytest.raises(KeyError, idx1.get_loc, 1.1)
        pytest.raises(TypeError, idx1.get_loc, idx1)

        # get the location of p1/p2 from
        # non-monotonic increasing/decreasing PeriodIndex with duplicate
        idx2 = pd.PeriodIndex([p2, p1, p2])
        expected_idx2_p1 = 1
        expected_idx2_p2 = np.array([True, False, True])

        assert idx2.get_loc(p1) == expected_idx2_p1
        assert idx2.get_loc(str(p1)) == expected_idx2_p1
        tm.assert_numpy_array_equal(idx2.get_loc(p2), expected_idx2_p2)
        tm.assert_numpy_array_equal(idx2.get_loc(str(p2)), expected_idx2_p2)

    def test_is_monotonic_increasing(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')

        idx_inc0 = pd.PeriodIndex([p0, p1, p2])
        idx_inc1 = pd.PeriodIndex([p0, p1, p1])
        idx_dec0 = pd.PeriodIndex([p2, p1, p0])
        idx_dec1 = pd.PeriodIndex([p2, p1, p1])
        idx = pd.PeriodIndex([p1, p2, p0])

        assert idx_inc0.is_monotonic_increasing
        assert idx_inc1.is_monotonic_increasing
        assert not idx_dec0.is_monotonic_increasing
        assert not idx_dec1.is_monotonic_increasing
        assert not idx.is_monotonic_increasing

    def test_is_monotonic_decreasing(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')

        idx_inc0 = pd.PeriodIndex([p0, p1, p2])
        idx_inc1 = pd.PeriodIndex([p0, p1, p1])
        idx_dec0 = pd.PeriodIndex([p2, p1, p0])
        idx_dec1 = pd.PeriodIndex([p2, p1, p1])
        idx = pd.PeriodIndex([p1, p2, p0])

        assert not idx_inc0.is_monotonic_decreasing
        assert not idx_inc1.is_monotonic_decreasing
        assert idx_dec0.is_monotonic_decreasing
        assert idx_dec1.is_monotonic_decreasing
        assert not idx.is_monotonic_decreasing

    def test_is_unique(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')

        idx0 = pd.PeriodIndex([p0, p1, p2])
        assert idx0.is_unique

        idx1 = pd.PeriodIndex([p1, p1, p2])
        assert not idx1.is_unique

    def test_contains(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')
        p3 = pd.Period('2017-09-04')

        ps0 = [p0, p1, p2]
        idx0 = pd.PeriodIndex(ps0)

        for p in ps0:
            assert idx0.contains(p)
            assert p in idx0

            assert idx0.contains(str(p))
            assert str(p) in idx0

        assert idx0.contains('2017-09-01 00:00:01')
        assert '2017-09-01 00:00:01' in idx0

        assert idx0.contains('2017-09')
        assert '2017-09' in idx0

        assert not idx0.contains(p3)
        assert p3 not in idx0

    def test_get_value(self):
        # GH 17717
        p0 = pd.Period('2017-09-01')
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')

        idx0 = pd.PeriodIndex([p0, p1, p2])
        input0 = np.array([1, 2, 3])
        expected0 = 2

        result0 = idx0.get_value(input0, p1)
        assert result0 == expected0

        idx1 = pd.PeriodIndex([p1, p1, p2])
        input1 = np.array([1, 2, 3])
        expected1 = np.array([1, 2])

        result1 = idx1.get_value(input1, p1)
        tm.assert_numpy_array_equal(result1, expected1)

        idx2 = pd.PeriodIndex([p1, p2, p1])
        input2 = np.array([1, 2, 3])
        expected2 = np.array([1, 3])

        result2 = idx2.get_value(input2, p1)
        tm.assert_numpy_array_equal(result2, expected2)

    def test_get_indexer(self):
        # GH 17717
        p1 = pd.Period('2017-09-01')
        p2 = pd.Period('2017-09-04')
        p3 = pd.Period('2017-09-07')

        tp0 = pd.Period('2017-08-31')
        tp1 = pd.Period('2017-09-02')
        tp2 = pd.Period('2017-09-05')
        tp3 = pd.Period('2017-09-09')

        idx = pd.PeriodIndex([p1, p2, p3])

        tm.assert_numpy_array_equal(idx.get_indexer(idx),
                                    np.array([0, 1, 2], dtype=np.intp))

        target = pd.PeriodIndex([tp0, tp1, tp2, tp3])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2, -1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 0, 1, 2], dtype=np.intp))

        res = idx.get_indexer(target, 'nearest',
                              tolerance=pd.Timedelta('1 day'))
        tm.assert_numpy_array_equal(res,
                                    np.array([0, 0, 1, -1], dtype=np.intp))

    def test_get_indexer_non_unique(self):
        # GH 17717
        p1 = pd.Period('2017-09-02')
        p2 = pd.Period('2017-09-03')
        p3 = pd.Period('2017-09-04')
        p4 = pd.Period('2017-09-05')

        idx1 = pd.PeriodIndex([p1, p2, p1])
        idx2 = pd.PeriodIndex([p2, p1, p3, p4])

        result = idx1.get_indexer_non_unique(idx2)
        expected_indexer = np.array([1, 0, 2, -1, -1], dtype=np.intp)
        expected_missing = np.array([2, 3], dtype=np.int64)

        tm.assert_numpy_array_equal(result[0], expected_indexer)
        tm.assert_numpy_array_equal(result[1], expected_missing)

    # TODO: This method came from test_period; de-dup with version above
    def test_get_loc2(self):
        idx = pd.period_range('2000-01-01', periods=3)

        for method in [None, 'pad', 'backfill', 'nearest']:
            assert idx.get_loc(idx[1], method) == 1
            assert idx.get_loc(idx[1].asfreq('H', how='start'), method) == 1
            assert idx.get_loc(idx[1].to_timestamp(), method) == 1
            assert idx.get_loc(idx[1].to_timestamp()
                               .to_pydatetime(), method) == 1
            assert idx.get_loc(str(idx[1]), method) == 1

        idx = pd.period_range('2000-01-01', periods=5)[::2]
        assert idx.get_loc('2000-01-02T12', method='nearest',
                           tolerance='1 day') == 1
        assert idx.get_loc('2000-01-02T12', method='nearest',
                           tolerance=pd.Timedelta('1D')) == 1
        assert idx.get_loc('2000-01-02T12', method='nearest',
                           tolerance=np.timedelta64(1, 'D')) == 1
        assert idx.get_loc('2000-01-02T12', method='nearest',
                           tolerance=timedelta(1)) == 1
        with tm.assert_raises_regex(ValueError,
                                    'unit abbreviation w/o a number'):
            idx.get_loc('2000-01-10', method='nearest', tolerance='foo')

        msg = 'Input has different freq from PeriodIndex\\(freq=D\\)'
        with tm.assert_raises_regex(ValueError, msg):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 hour')
        with pytest.raises(KeyError):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 day')
        with pytest.raises(
                ValueError,
                match='list-like tolerance size must match target index size'):
            idx.get_loc('2000-01-10', method='nearest',
                        tolerance=[pd.Timedelta('1 day').to_timedelta64(),
                                   pd.Timedelta('1 day').to_timedelta64()])

    # TODO: This method came from test_period; de-dup with version above
    def test_get_indexer2(self):
        idx = pd.period_range('2000-01-01', periods=3).asfreq('H', how='start')
        tm.assert_numpy_array_equal(idx.get_indexer(idx),
                                    np.array([0, 1, 2], dtype=np.intp))

        target = pd.PeriodIndex(['1999-12-31T23', '2000-01-01T12',
                                 '2000-01-02T01'], freq='H')
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 1, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest',
                                                    tolerance='1 hour'),
                                    np.array([0, -1, 1], dtype=np.intp))

        msg = 'Input has different freq from PeriodIndex\\(freq=H\\)'
        with tm.assert_raises_regex(ValueError, msg):
            idx.get_indexer(target, 'nearest', tolerance='1 minute')

        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest',
                                                    tolerance='1 day'),
                                    np.array([0, 1, 1], dtype=np.intp))
        tol_raw = [pd.Timedelta('1 hour'),
                   pd.Timedelta('1 hour'),
                   np.timedelta64(1, 'D'), ]
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest',
                            tolerance=[np.timedelta64(x) for x in tol_raw]),
            np.array([0, -1, 1], dtype=np.intp))
        tol_bad = [pd.Timedelta('2 hour').to_timedelta64(),
                   pd.Timedelta('1 hour').to_timedelta64(),
                   np.timedelta64(1, 'M'), ]
        with pytest.raises(
                libperiod.IncompatibleFrequency,
                match='Input has different freq from'):
            idx.get_indexer(target, 'nearest', tolerance=tol_bad)

    def test_indexing(self):
        # GH 4390, iat incorrectly indexing
        index = period_range('1/1/2001', periods=10)
        s = Series(np.random.randn(10), index=index)
        expected = s[index[0]]
        result = s.iat[0]
        assert expected == result

    def test_period_index_indexer(self):
        # GH4125
        idx = pd.period_range('2002-01', '2003-12', freq='M')
        df = pd.DataFrame(pd.np.random.randn(24, 10), index=idx)
        tm.assert_frame_equal(df, df.loc[idx])
        tm.assert_frame_equal(df, df.loc[list(idx)])
        tm.assert_frame_equal(df, df.loc[list(idx)])
        tm.assert_frame_equal(df.iloc[0:5], df.loc[idx[0:5]])
        tm.assert_frame_equal(df, df.loc[list(idx)])
