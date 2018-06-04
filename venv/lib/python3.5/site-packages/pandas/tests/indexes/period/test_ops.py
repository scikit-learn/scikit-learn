
import numpy as np
import pytest

import pandas as pd
import pandas._libs.tslib as tslib
import pandas.util.testing as tm
from pandas import (DatetimeIndex, PeriodIndex, Series, Period,
                    _np_version_under1p10, Index)

from pandas.tests.test_base import Ops


class TestPeriodIndexOps(Ops):

    def setup_method(self, method):
        super(TestPeriodIndexOps, self).setup_method(method)
        mask = lambda x: (isinstance(x, DatetimeIndex) or
                          isinstance(x, PeriodIndex))
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = [o for o in self.objs if not mask(o)]

    def test_ops_properties(self):
        f = lambda x: isinstance(x, PeriodIndex)
        self.check_ops_properties(PeriodIndex._field_ops, f)
        self.check_ops_properties(PeriodIndex._object_ops, f)
        self.check_ops_properties(PeriodIndex._bool_ops, f)

    def test_minmax(self):

        # monotonic
        idx1 = pd.PeriodIndex([pd.NaT, '2011-01-01', '2011-01-02',
                               '2011-01-03'], freq='D')
        assert idx1.is_monotonic

        # non-monotonic
        idx2 = pd.PeriodIndex(['2011-01-01', pd.NaT, '2011-01-03',
                               '2011-01-02', pd.NaT], freq='D')
        assert not idx2.is_monotonic

        for idx in [idx1, idx2]:
            assert idx.min() == pd.Period('2011-01-01', freq='D')
            assert idx.max() == pd.Period('2011-01-03', freq='D')
        assert idx1.argmin() == 1
        assert idx2.argmin() == 0
        assert idx1.argmax() == 3
        assert idx2.argmax() == 2

        for op in ['min', 'max']:
            # Return NaT
            obj = PeriodIndex([], freq='M')
            result = getattr(obj, op)()
            assert result is tslib.NaT

            obj = PeriodIndex([pd.NaT], freq='M')
            result = getattr(obj, op)()
            assert result is tslib.NaT

            obj = PeriodIndex([pd.NaT, pd.NaT, pd.NaT], freq='M')
            result = getattr(obj, op)()
            assert result is tslib.NaT

    def test_numpy_minmax(self):
        pr = pd.period_range(start='2016-01-15', end='2016-01-20')

        assert np.min(pr) == Period('2016-01-15', freq='D')
        assert np.max(pr) == Period('2016-01-20', freq='D')

        errmsg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, errmsg, np.min, pr, out=0)
        tm.assert_raises_regex(ValueError, errmsg, np.max, pr, out=0)

        assert np.argmin(pr) == 0
        assert np.argmax(pr) == 5

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmin, pr, out=0)
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmax, pr, out=0)

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H',
                                   'T', 'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day',
                                   'hour', 'minute', 'second',
                                   'millisecond', 'microsecond']):

            idx = pd.period_range(start='2013-04-01', periods=30, freq=freq)
            assert idx.resolution == expected

    def test_value_counts_unique(self):
        # GH 7735
        idx = pd.period_range('2011-01-01 09:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = PeriodIndex(np.repeat(idx.values, range(1, len(idx) + 1)),
                          freq='H')

        exp_idx = PeriodIndex(['2011-01-01 18:00', '2011-01-01 17:00',
                               '2011-01-01 16:00', '2011-01-01 15:00',
                               '2011-01-01 14:00', '2011-01-01 13:00',
                               '2011-01-01 12:00', '2011-01-01 11:00',
                               '2011-01-01 10:00',
                               '2011-01-01 09:00'], freq='H')
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = pd.period_range('2011-01-01 09:00', freq='H',
                                   periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 09:00',
                           '2013-01-01 09:00', '2013-01-01 08:00',
                           '2013-01-01 08:00', pd.NaT], freq='H')

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00'],
                              freq='H')
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00',
                               pd.NaT], freq='H')
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.period_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        result = idx.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

        idx_dup = idx.append(idx)  # freq will not be reset
        result = idx_dup.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')
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

    def test_order_compat(self):
        def _check_freq(index, expected_index):
            if isinstance(index, PeriodIndex):
                assert index.freq == expected_index.freq

        pidx = PeriodIndex(['2011', '2012', '2013'], name='pidx', freq='A')
        # for compatibility check
        iidx = Index([2011, 2012, 2013], name='idx')
        for idx in [pidx, iidx]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, idx)
            _check_freq(ordered, idx)

            ordered = idx.sort_values(ascending=False)
            tm.assert_index_equal(ordered, idx[::-1])
            _check_freq(ordered, idx[::-1])

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, idx)
            tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]),
                                        check_dtype=False)
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            tm.assert_index_equal(ordered, idx[::-1])
            tm.assert_numpy_array_equal(indexer, np.array([2, 1, 0]),
                                        check_dtype=False)
            _check_freq(ordered, idx[::-1])

        pidx = PeriodIndex(['2011', '2013', '2015', '2012',
                            '2011'], name='pidx', freq='A')
        pexpected = PeriodIndex(
            ['2011', '2011', '2012', '2013', '2015'], name='pidx', freq='A')
        # for compatibility check
        iidx = Index([2011, 2013, 2015, 2012, 2011], name='idx')
        iexpected = Index([2011, 2011, 2012, 2013, 2015], name='idx')
        for idx, expected in [(pidx, pexpected), (iidx, iexpected)]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, expected)
            _check_freq(ordered, idx)

            ordered = idx.sort_values(ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])

            exp = np.array([2, 1, 3, 4, 0])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            _check_freq(ordered, idx)

        pidx = PeriodIndex(['2011', '2013', 'NaT', '2011'], name='pidx',
                           freq='D')

        result = pidx.sort_values()
        expected = PeriodIndex(['NaT', '2011', '2011', '2013'],
                               name='pidx', freq='D')
        tm.assert_index_equal(result, expected)
        assert result.freq == 'D'

        result = pidx.sort_values(ascending=False)
        expected = PeriodIndex(
            ['2013', '2011', '2011', 'NaT'], name='pidx', freq='D')
        tm.assert_index_equal(result, expected)
        assert result.freq == 'D'

    def test_order(self):
        for freq in ['D', '2D', '4D']:
            idx = PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03'],
                              freq=freq, name='idx')

            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, idx)
            assert ordered.freq == idx.freq

            ordered = idx.sort_values(ascending=False)
            expected = idx[::-1]
            tm.assert_index_equal(ordered, expected)
            assert ordered.freq == expected.freq
            assert ordered.freq == freq

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, idx)
            tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]),
                                        check_dtype=False)
            assert ordered.freq == idx.freq
            assert ordered.freq == freq

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            expected = idx[::-1]
            tm.assert_index_equal(ordered, expected)
            tm.assert_numpy_array_equal(indexer, np.array([2, 1, 0]),
                                        check_dtype=False)
            assert ordered.freq == expected.freq
            assert ordered.freq == freq

        idx1 = PeriodIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                            '2011-01-02', '2011-01-01'], freq='D', name='idx1')
        exp1 = PeriodIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                            '2011-01-03', '2011-01-05'], freq='D', name='idx1')

        idx2 = PeriodIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                            '2011-01-02', '2011-01-01'],
                           freq='D', name='idx2')
        exp2 = PeriodIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                            '2011-01-03', '2011-01-05'],
                           freq='D', name='idx2')

        idx3 = PeriodIndex([pd.NaT, '2011-01-03', '2011-01-05',
                            '2011-01-02', pd.NaT], freq='D', name='idx3')
        exp3 = PeriodIndex([pd.NaT, pd.NaT, '2011-01-02', '2011-01-03',
                            '2011-01-05'], freq='D', name='idx3')

        for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3)]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, expected)
            assert ordered.freq == 'D'

            ordered = idx.sort_values(ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])
            assert ordered.freq == 'D'

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            assert ordered.freq == 'D'

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])

            exp = np.array([2, 1, 3, 4, 0])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            assert ordered.freq == 'D'

    def test_nat_new(self):

        idx = pd.period_range('2011-01', freq='M', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.PeriodIndex([pd.NaT] * 5, freq='M', name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # This is tested in test_arithmetic
        pass

    def test_repeat(self):
        index = pd.period_range('2001-01-01', periods=2, freq='D')
        exp = pd.PeriodIndex(['2001-01-01', '2001-01-01',
                              '2001-01-02', '2001-01-02'], freq='D')
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)

        index = pd.period_range('2001-01-01', periods=2, freq='2D')
        exp = pd.PeriodIndex(['2001-01-01', '2001-01-01',
                              '2001-01-03', '2001-01-03'], freq='2D')
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)

        index = pd.PeriodIndex(['2001-01', 'NaT', '2003-01'], freq='M')
        exp = pd.PeriodIndex(['2001-01', '2001-01', '2001-01',
                              'NaT', 'NaT', 'NaT',
                              '2003-01', '2003-01', '2003-01'], freq='M')
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)

    def test_nat(self):
        assert pd.PeriodIndex._na_value is pd.NaT
        assert pd.PeriodIndex([], freq='M')._na_value is pd.NaT

        idx = pd.PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        assert not idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([], dtype=np.intp))

        idx = pd.PeriodIndex(['2011-01-01', 'NaT'], freq='D')
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        assert idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([1], dtype=np.intp))

    @pytest.mark.parametrize('freq', ['D', 'M'])
    def test_equals(self, freq):
        # GH#13107
        idx = pd.PeriodIndex(['2011-01-01', '2011-01-02', 'NaT'],
                             freq=freq)
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))
        assert idx.astype(object).equals(idx)
        assert idx.astype(object).equals(idx.astype(object))
        assert not idx.equals(list(idx))
        assert not idx.equals(pd.Series(idx))

        idx2 = pd.PeriodIndex(['2011-01-01', '2011-01-02', 'NaT'],
                              freq='H')
        assert not idx.equals(idx2)
        assert not idx.equals(idx2.copy())
        assert not idx.equals(idx2.astype(object))
        assert not idx.astype(object).equals(idx2)
        assert not idx.equals(list(idx2))
        assert not idx.equals(pd.Series(idx2))

        # same internal, different tz
        idx3 = pd.PeriodIndex._simple_new(idx.asi8, freq='H')
        tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
        assert not idx.equals(idx3)
        assert not idx.equals(idx3.copy())
        assert not idx.equals(idx3.astype(object))
        assert not idx.astype(object).equals(idx3)
        assert not idx.equals(list(idx3))
        assert not idx.equals(pd.Series(idx3))

    def test_freq_setter_deprecated(self):
        # GH 20678
        idx = pd.period_range('2018Q1', periods=4, freq='Q')

        # no warning for getter
        with tm.assert_produces_warning(None):
            idx.freq

        # warning for setter
        with tm.assert_produces_warning(FutureWarning):
            idx.freq = pd.offsets.Day()


class TestPeriodIndexSeriesMethods(object):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        if isinstance(expected, pd.Index):
            tm.assert_index_equal(result, expected)
        else:
            # comp op results in bool
            tm.assert_numpy_array_equal(result, expected)

        s = pd.Series(values)
        result = func(s)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_comp_period(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        f = lambda x: x == pd.Period('2011-03', freq='M')
        exp = np.array([False, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') == x
        self._check(idx, f, exp)

        f = lambda x: x != pd.Period('2011-03', freq='M')
        exp = np.array([True, True, False, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') != x
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, True, True, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x > pd.Period('2011-03', freq='M')
        exp = np.array([False, False, False, True], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, True, True, False], dtype=np.bool)
        self._check(idx, f, exp)

    def test_pi_comp_period_nat(self):
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03',
                           '2011-04'], freq='M', name='idx')

        f = lambda x: x == pd.Period('2011-03', freq='M')
        exp = np.array([False, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') == x
        self._check(idx, f, exp)

        f = lambda x: x == tslib.NaT
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: tslib.NaT == x
        self._check(idx, f, exp)

        f = lambda x: x != pd.Period('2011-03', freq='M')
        exp = np.array([True, True, False, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') != x
        self._check(idx, f, exp)

        f = lambda x: x != tslib.NaT
        exp = np.array([True, True, True, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: tslib.NaT != x
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x < pd.Period('2011-03', freq='M')
        exp = np.array([True, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x > tslib.NaT
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: tslib.NaT >= x
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)
