import pytest
import warnings
import numpy as np
from datetime import datetime

import pandas as pd
import pandas._libs.tslib as tslib
import pandas.util.testing as tm
from pandas import (DatetimeIndex, PeriodIndex, Series, Timestamp,
                    date_range, _np_version_under1p10, Index,
                    bdate_range)
from pandas.tseries.offsets import BMonthEnd, CDay, BDay, Day, Hour
from pandas.tests.test_base import Ops
from pandas.core.dtypes.generic import ABCDateOffset


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
                        'dateutil/Asia/Singapore',
                        'dateutil/US/Pacific'])
def tz_fixture(request):
    return request.param


START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)


class TestDatetimeIndexOps(Ops):

    def setup_method(self, method):
        super(TestDatetimeIndexOps, self).setup_method(method)
        mask = lambda x: (isinstance(x, DatetimeIndex) or
                          isinstance(x, PeriodIndex))
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = [o for o in self.objs if not mask(o)]

    def test_ops_properties(self):
        f = lambda x: isinstance(x, DatetimeIndex)
        self.check_ops_properties(DatetimeIndex._field_ops, f)
        self.check_ops_properties(DatetimeIndex._object_ops, f)
        self.check_ops_properties(DatetimeIndex._bool_ops, f)

    def test_ops_properties_basic(self):

        # sanity check that the behavior didn't change
        # GH7206
        for op in ['year', 'day', 'second', 'weekday']:
            pytest.raises(TypeError, lambda x: getattr(self.dt_series, op))

        # attribute access should still work!
        s = Series(dict(year=2000, month=1, day=10))
        assert s.year == 2000
        assert s.month == 1
        assert s.day == 10
        pytest.raises(AttributeError, lambda: s.weekday)

    def test_minmax_tz(self, tz_fixture):
        tz = tz_fixture
        # monotonic
        idx1 = pd.DatetimeIndex(['2011-01-01', '2011-01-02',
                                 '2011-01-03'], tz=tz)
        assert idx1.is_monotonic

        # non-monotonic
        idx2 = pd.DatetimeIndex(['2011-01-01', pd.NaT, '2011-01-03',
                                 '2011-01-02', pd.NaT], tz=tz)
        assert not idx2.is_monotonic

        for idx in [idx1, idx2]:
            assert idx.min() == Timestamp('2011-01-01', tz=tz)
            assert idx.max() == Timestamp('2011-01-03', tz=tz)
            assert idx.argmin() == 0
            assert idx.argmax() == 2

    @pytest.mark.parametrize('op', ['min', 'max'])
    def test_minmax_nat(self, op):
        # Return NaT
        obj = DatetimeIndex([])
        assert pd.isna(getattr(obj, op)())

        obj = DatetimeIndex([pd.NaT])
        assert pd.isna(getattr(obj, op)())

        obj = DatetimeIndex([pd.NaT, pd.NaT, pd.NaT])
        assert pd.isna(getattr(obj, op)())

    def test_numpy_minmax(self):
        dr = pd.date_range(start='2016-01-15', end='2016-01-20')

        assert np.min(dr) == Timestamp('2016-01-15 00:00:00', freq='D')
        assert np.max(dr) == Timestamp('2016-01-20 00:00:00', freq='D')

        errmsg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, errmsg, np.min, dr, out=0)
        tm.assert_raises_regex(ValueError, errmsg, np.max, dr, out=0)

        assert np.argmin(dr) == 0
        assert np.argmax(dr) == 5

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmin, dr, out=0)
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmax, dr, out=0)

    def test_repeat_range(self, tz_fixture):
        tz = tz_fixture
        rng = date_range('1/1/2000', '1/1/2001')

        result = rng.repeat(5)
        assert result.freq is None
        assert len(result) == 5 * len(rng)

        index = pd.date_range('2001-01-01', periods=2, freq='D', tz=tz)
        exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                '2001-01-02', '2001-01-02'], tz=tz)
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = pd.date_range('2001-01-01', periods=2, freq='2D', tz=tz)
        exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                '2001-01-03', '2001-01-03'], tz=tz)
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = pd.DatetimeIndex(['2001-01-01', 'NaT', '2003-01-01'],
                                 tz=tz)
        exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01', '2001-01-01',
                                'NaT', 'NaT', 'NaT',
                                '2003-01-01', '2003-01-01', '2003-01-01'],
                               tz=tz)
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

    def test_repeat(self, tz_fixture):
        tz = tz_fixture
        reps = 2
        msg = "the 'axis' parameter is not supported"

        rng = pd.date_range(start='2016-01-01', periods=2,
                            freq='30Min', tz=tz)

        expected_rng = DatetimeIndex([
            Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
            Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
            Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
            Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
        ])

        res = rng.repeat(reps)
        tm.assert_index_equal(res, expected_rng)
        assert res.freq is None

        tm.assert_index_equal(np.repeat(rng, reps), expected_rng)
        tm.assert_raises_regex(ValueError, msg, np.repeat,
                               rng, reps, axis=1)

    def test_resolution(self, tz_fixture):
        tz = tz_fixture
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T',
                                   'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day', 'hour',
                                   'minute', 'second', 'millisecond',
                                   'microsecond']):
            idx = pd.date_range(start='2013-04-01', periods=30, freq=freq,
                                tz=tz)
            assert idx.resolution == expected

    def test_value_counts_unique(self, tz_fixture):
        tz = tz_fixture
        # GH 7735
        idx = pd.date_range('2011-01-01 09:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = DatetimeIndex(np.repeat(idx.values, range(1, len(idx) + 1)),
                            tz=tz)

        exp_idx = pd.date_range('2011-01-01 18:00', freq='-1H', periods=10,
                                tz=tz)
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = pd.date_range('2011-01-01 09:00', freq='H', periods=10,
                                 tz=tz)
        tm.assert_index_equal(idx.unique(), expected)

        idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 09:00',
                             '2013-01-01 09:00', '2013-01-01 08:00',
                             '2013-01-01 08:00', pd.NaT], tz=tz)

        exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00'],
                                tz=tz)
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00',
                                 pd.NaT], tz=tz)
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False),
                                   expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(DatetimeIndex,
                       ([0, 1, 0], [0, 0, -1], [0, -1, -1],
                        ['2015', '2015', '2016'], ['2015', '2015', '2014'])):
            assert idx[0] in idx

    @pytest.mark.parametrize('idx',
                             [
                                 DatetimeIndex(
                                     ['2011-01-01',
                                      '2011-01-02',
                                      '2011-01-03'],
                                     freq='D', name='idx'),
                                 DatetimeIndex(
                                     ['2011-01-01 09:00',
                                      '2011-01-01 10:00',
                                      '2011-01-01 11:00'],
                                     freq='H', name='tzidx', tz='Asia/Tokyo')
                             ])
    def test_order_with_freq(self, idx):
        ordered = idx.sort_values()
        tm.assert_index_equal(ordered, idx)
        assert ordered.freq == idx.freq

        ordered = idx.sort_values(ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        assert ordered.freq == expected.freq
        assert ordered.freq.n == -1

        ordered, indexer = idx.sort_values(return_indexer=True)
        tm.assert_index_equal(ordered, idx)
        tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]),
                                    check_dtype=False)
        assert ordered.freq == idx.freq

        ordered, indexer = idx.sort_values(return_indexer=True,
                                           ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        tm.assert_numpy_array_equal(indexer,
                                    np.array([2, 1, 0]),
                                    check_dtype=False)
        assert ordered.freq == expected.freq
        assert ordered.freq.n == -1

    @pytest.mark.parametrize('index_dates,expected_dates', [
        (['2011-01-01', '2011-01-03', '2011-01-05',
          '2011-01-02', '2011-01-01'],
         ['2011-01-01', '2011-01-01', '2011-01-02',
          '2011-01-03', '2011-01-05']),
        (['2011-01-01', '2011-01-03', '2011-01-05',
          '2011-01-02', '2011-01-01'],
         ['2011-01-01', '2011-01-01', '2011-01-02',
          '2011-01-03', '2011-01-05']),
        ([pd.NaT, '2011-01-03', '2011-01-05',
          '2011-01-02', pd.NaT],
         [pd.NaT, pd.NaT, '2011-01-02', '2011-01-03',
          '2011-01-05'])
    ])
    def test_order_without_freq(self, index_dates, expected_dates, tz_fixture):
        tz = tz_fixture

        # without freq
        index = DatetimeIndex(index_dates, tz=tz, name='idx')
        expected = DatetimeIndex(expected_dates, tz=tz, name='idx')

        ordered = index.sort_values()
        tm.assert_index_equal(ordered, expected)
        assert ordered.freq is None

        ordered = index.sort_values(ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])
        assert ordered.freq is None

        ordered, indexer = index.sort_values(return_indexer=True)
        tm.assert_index_equal(ordered, expected)

        exp = np.array([0, 4, 3, 1, 2])
        tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
        assert ordered.freq is None

        ordered, indexer = index.sort_values(return_indexer=True,
                                             ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])

        exp = np.array([2, 1, 3, 4, 0])
        tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
        assert ordered.freq is None

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        result = idx.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

        idx_dup = idx.append(idx)
        assert idx_dup.freq is None  # freq is reset
        result = idx_dup.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert result.freq is None

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        idx = base.append(base[:5])

        res = idx.drop_duplicates()
        tm.assert_index_equal(res, base)
        res = Series(idx).drop_duplicates()
        tm.assert_series_equal(res, Series(base))

        res = idx.drop_duplicates(keep='last')
        exp = base[5:].append(base[:5])
        tm.assert_index_equal(res, exp)
        res = Series(idx).drop_duplicates(keep='last')
        tm.assert_series_equal(res, Series(exp, index=np.arange(5, 36)))

        res = idx.drop_duplicates(keep=False)
        tm.assert_index_equal(res, base[5:])
        res = Series(idx).drop_duplicates(keep=False)
        tm.assert_series_equal(res, Series(base[5:], index=np.arange(5, 31)))

    @pytest.mark.parametrize('freq', [
        'A', '2A', '-2A', 'Q', '-1Q', 'M', '-1M', 'D', '3D',
        '-3D', 'W', '-1W', 'H', '2H', '-2H', 'T', '2T', 'S',
        '-3S'])
    def test_infer_freq(self, freq):
        # GH 11018
        idx = pd.date_range('2011-01-01 09:00:00', freq=freq, periods=10)
        result = pd.DatetimeIndex(idx.asi8, freq='infer')
        tm.assert_index_equal(idx, result)
        assert result.freq == freq

    def test_nat_new(self):
        idx = pd.date_range('2011-01-01', freq='D', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.DatetimeIndex([pd.NaT] * 5, name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_nat(self, tz_naive_fixture):
        timezone = tz_naive_fixture
        assert pd.DatetimeIndex._na_value is pd.NaT
        assert pd.DatetimeIndex([])._na_value is pd.NaT

        idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], tz=timezone)
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        assert not idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([], dtype=np.intp))

        idx = pd.DatetimeIndex(['2011-01-01', 'NaT'], tz=timezone)
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        assert idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'])
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))
        assert idx.astype(object).equals(idx)
        assert idx.astype(object).equals(idx.astype(object))
        assert not idx.equals(list(idx))
        assert not idx.equals(pd.Series(idx))

        idx2 = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'],
                                tz='US/Pacific')
        assert not idx.equals(idx2)
        assert not idx.equals(idx2.copy())
        assert not idx.equals(idx2.astype(object))
        assert not idx.astype(object).equals(idx2)
        assert not idx.equals(list(idx2))
        assert not idx.equals(pd.Series(idx2))

        # same internal, different tz
        idx3 = pd.DatetimeIndex._simple_new(idx.asi8, tz='US/Pacific')
        tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
        assert not idx.equals(idx3)
        assert not idx.equals(idx3.copy())
        assert not idx.equals(idx3.astype(object))
        assert not idx.astype(object).equals(idx3)
        assert not idx.equals(list(idx3))
        assert not idx.equals(pd.Series(idx3))

    @pytest.mark.parametrize('values', [
        ['20180101', '20180103', '20180105'], []])
    @pytest.mark.parametrize('freq', [
        '2D', Day(2), '2B', BDay(2), '48H', Hour(48)])
    @pytest.mark.parametrize('tz', [None, 'US/Eastern'])
    def test_freq_setter(self, values, freq, tz):
        # GH 20678
        idx = DatetimeIndex(values, tz=tz)

        # can set to an offset, converting from string if necessary
        idx.freq = freq
        assert idx.freq == freq
        assert isinstance(idx.freq, ABCDateOffset)

        # can reset to None
        idx.freq = None
        assert idx.freq is None

    def test_freq_setter_errors(self):
        # GH 20678
        idx = DatetimeIndex(['20180101', '20180103', '20180105'])

        # setting with an incompatible freq
        msg = ('Inferred frequency 2D from passed values does not conform to '
               'passed frequency 5D')
        with tm.assert_raises_regex(ValueError, msg):
            idx.freq = '5D'

        # setting with non-freq string
        with tm.assert_raises_regex(ValueError, 'Invalid frequency'):
            idx.freq = 'foo'

    def test_offset_deprecated(self):
        # GH 20716
        idx = pd.DatetimeIndex(['20180101', '20180102'])

        # getter deprecated
        with tm.assert_produces_warning(FutureWarning):
            idx.offset

        # setter deprecated
        with tm.assert_produces_warning(FutureWarning):
            idx.offset = BDay()


class TestBusinessDatetimeIndex(object):

    def setup_method(self, method):
        self.rng = bdate_range(START, END)

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_pickle_unpickle(self):
        unpickled = tm.round_trip_pickle(self.rng)
        assert unpickled.freq is not None

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_shift(self):
        shifted = self.rng.shift(5)
        assert shifted[0] == self.rng[5]
        assert shifted.freq == self.rng.freq

        shifted = self.rng.shift(-5)
        assert shifted[5] == self.rng[0]
        assert shifted.freq == self.rng.freq

        shifted = self.rng.shift(0)
        assert shifted[0] == self.rng[0]
        assert shifted.freq == self.rng.freq

        rng = date_range(START, END, freq=BMonthEnd())
        shifted = rng.shift(1, freq=BDay())
        assert shifted[0] == rng[0] + BDay()

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))

    def test_identical(self):
        t1 = self.rng.copy()
        t2 = self.rng.copy()
        assert t1.identical(t2)

        # name
        t1 = t1.rename('foo')
        assert t1.equals(t2)
        assert not t1.identical(t2)
        t2 = t2.rename('foo')
        assert t1.identical(t2)

        # freq
        t2v = Index(t2.values)
        assert t1.equals(t2v)
        assert not t1.identical(t2v)


class TestCustomDatetimeIndex(object):
    def setup_method(self, method):
        self.rng = bdate_range(START, END, freq='C')

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_shift(self):

        shifted = self.rng.shift(5)
        assert shifted[0] == self.rng[5]
        assert shifted.freq == self.rng.freq

        shifted = self.rng.shift(-5)
        assert shifted[5] == self.rng[0]
        assert shifted.freq == self.rng.freq

        shifted = self.rng.shift(0)
        assert shifted[0] == self.rng[0]
        assert shifted.freq == self.rng.freq

        # PerformanceWarning
        with warnings.catch_warnings(record=True):
            rng = date_range(START, END, freq=BMonthEnd())
            shifted = rng.shift(1, freq=CDay())
            assert shifted[0] == rng[0] + CDay()

    def test_pickle_unpickle(self):
        unpickled = tm.round_trip_pickle(self.rng)
        assert unpickled.freq is not None

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))
