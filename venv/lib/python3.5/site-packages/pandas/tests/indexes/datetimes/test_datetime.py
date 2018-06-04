import warnings

import pytest

import numpy as np
from datetime import date

import dateutil
import pandas as pd
import pandas.util.testing as tm
from pandas.compat import lrange
from pandas import (DatetimeIndex, Index, date_range, DataFrame,
                    Timestamp, offsets)

from pandas.util.testing import assert_almost_equal

randn = np.random.randn


class TestDatetimeIndex(object):

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = date_range('20130101', periods=3, tz='US/Eastern', name='foo')
        unpickled = tm.round_trip_pickle(index)
        tm.assert_index_equal(index, unpickled)

    def test_reindex_preserves_tz_if_target_is_empty_list_or_array(self):
        # GH7774
        index = date_range('20130101', periods=3, tz='US/Eastern')
        assert str(index.reindex([])[0].tz) == 'US/Eastern'
        assert str(index.reindex(np.array([]))[0].tz) == 'US/Eastern'

    def test_time_loc(self):  # GH8667
        from datetime import time
        from pandas._libs.index import _SIZE_CUTOFF

        ns = _SIZE_CUTOFF + np.array([-100, 100], dtype=np.int64)
        key = time(15, 11, 30)
        start = key.hour * 3600 + key.minute * 60 + key.second
        step = 24 * 3600

        for n in ns:
            idx = pd.date_range('2014-11-26', periods=n, freq='S')
            ts = pd.Series(np.random.randn(n), index=idx)
            i = np.arange(start, n, step)

            tm.assert_numpy_array_equal(ts.index.get_loc(key), i,
                                        check_dtype=False)
            tm.assert_series_equal(ts[key], ts.iloc[i])

            left, right = ts.copy(), ts.copy()
            left[key] *= -10
            right.iloc[i] *= -10
            tm.assert_series_equal(left, right)

    def test_time_overflow_for_32bit_machines(self):
        # GH8943.  On some machines NumPy defaults to np.int32 (for example,
        # 32-bit Linux machines).  In the function _generate_regular_range
        # found in tseries/index.py, `periods` gets multiplied by `strides`
        # (which has value 1e9) and since the max value for np.int32 is ~2e9,
        # and since those machines won't promote np.int32 to np.int64, we get
        # overflow.
        periods = np.int_(1000)

        idx1 = pd.date_range(start='2000', periods=periods, freq='S')
        assert len(idx1) == periods

        idx2 = pd.date_range(end='2000', periods=periods, freq='S')
        assert len(idx2) == periods

    def test_nat(self):
        assert DatetimeIndex([np.nan])[0] is pd.NaT

    def test_week_of_month_frequency(self):
        # GH 5348: "ValueError: Could not evaluate WOM-1SUN" shouldn't raise
        d1 = date(2002, 9, 1)
        d2 = date(2013, 10, 27)
        d3 = date(2012, 9, 30)
        idx1 = DatetimeIndex([d1, d2])
        idx2 = DatetimeIndex([d3])
        result_append = idx1.append(idx2)
        expected = DatetimeIndex([d1, d2, d3])
        tm.assert_index_equal(result_append, expected)
        result_union = idx1.union(idx2)
        expected = DatetimeIndex([d1, d3, d2])
        tm.assert_index_equal(result_union, expected)

        # GH 5115
        result = date_range("2013-1-1", periods=4, freq='WOM-1SAT')
        dates = ['2013-01-05', '2013-02-02', '2013-03-02', '2013-04-06']
        expected = DatetimeIndex(dates, freq='WOM-1SAT')
        tm.assert_index_equal(result, expected)

    def test_hash_error(self):
        index = date_range('20010101', periods=10)
        with tm.assert_raises_regex(TypeError, "unhashable type: %r" %
                                    type(index).__name__):
            hash(index)

    def test_stringified_slice_with_tz(self):
        # GH2658
        import datetime
        start = datetime.datetime.now()
        idx = DatetimeIndex(start=start, freq="1d", periods=10)
        df = DataFrame(lrange(10), index=idx)
        df["2013-01-14 23:44:34.437768-05:00":]  # no exception here

    def test_append_join_nondatetimeindex(self):
        rng = date_range('1/1/2000', periods=10)
        idx = Index(['a', 'b', 'c', 'd'])

        result = rng.append(idx)
        assert isinstance(result[0], Timestamp)

        # it works
        rng.join(idx, how='outer')

    def test_map(self):
        rng = date_range('1/1/2000', periods=10)

        f = lambda x: x.strftime('%Y%m%d')
        result = rng.map(f)
        exp = Index([f(x) for x in rng], dtype='<U8')
        tm.assert_index_equal(result, exp)

    def test_iteration_preserves_tz(self):
        # see gh-8890
        index = date_range("2012-01-01", periods=3, freq='H', tz='US/Eastern')

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result == expected

        index = date_range("2012-01-01", periods=3, freq='H',
                           tz=dateutil.tz.tzoffset(None, -28800))

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result._repr_base == expected._repr_base
            assert result == expected

        # 9100
        index = pd.DatetimeIndex(['2014-12-01 03:32:39.987000-08:00',
                                  '2014-12-01 04:12:34.987000-08:00'])
        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            assert result._repr_base == expected._repr_base
            assert result == expected

    @pytest.mark.parametrize('periods', [0, 9999, 10000, 10001])
    def test_iteration_over_chunksize(self, periods):
        # GH21012

        index = date_range('2000-01-01 00:00:00', periods=periods, freq='min')
        num = 0
        for stamp in index:
            assert index[num] == stamp
            num += 1
        assert num == len(index)

    def test_misc_coverage(self):
        rng = date_range('1/1/2000', periods=5)
        result = rng.groupby(rng.day)
        assert isinstance(list(result.values())[0][0], Timestamp)

        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        assert not idx.equals(list(idx))

        non_datetime = Index(list('abc'))
        assert not idx.equals(list(non_datetime))

    def test_string_index_series_name_converted(self):
        # #1644
        df = DataFrame(np.random.randn(10, 4),
                       index=date_range('1/1/2000', periods=10))

        result = df.loc['1/3/2000']
        assert result.name == df.index[2]

        result = df.T['1/3/2000']
        assert result.name == df.index[2]

    def test_get_duplicates(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-02',
                             '2000-01-03', '2000-01-03', '2000-01-04'])

        with warnings.catch_warnings(record=True):
            # Deprecated - see GH20239
            result = idx.get_duplicates()

        ex = DatetimeIndex(['2000-01-02', '2000-01-03'])
        tm.assert_index_equal(result, ex)

    def test_argmin_argmax(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])
        assert idx.argmin() == 1
        assert idx.argmax() == 0

    def test_sort_values(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])

        ordered = idx.sort_values()
        assert ordered.is_monotonic

        ordered = idx.sort_values(ascending=False)
        assert ordered[::-1].is_monotonic

        ordered, dexer = idx.sort_values(return_indexer=True)
        assert ordered.is_monotonic
        tm.assert_numpy_array_equal(dexer, np.array([1, 2, 0], dtype=np.intp))

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        assert ordered[::-1].is_monotonic
        tm.assert_numpy_array_equal(dexer, np.array([0, 2, 1], dtype=np.intp))

    def test_map_bug_1677(self):
        index = DatetimeIndex(['2012-04-25 09:30:00.393000'])
        f = index.asof

        result = index.map(f)
        expected = Index([f(index[0])])
        tm.assert_index_equal(result, expected)

    def test_groupby_function_tuple_1677(self):
        df = DataFrame(np.random.rand(100),
                       index=date_range("1/1/2000", periods=100))
        monthly_group = df.groupby(lambda x: (x.year, x.month))

        result = monthly_group.mean()
        assert isinstance(result.index[0], tuple)

    def test_append_numpy_bug_1681(self):
        # another datetime64 bug
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        a = DataFrame()
        c = DataFrame({'A': 'foo', 'B': dr}, index=dr)

        result = a.append(c)
        assert (result['B'] == dr).all()

    def test_isin(self):
        index = tm.makeDateIndex(4)
        result = index.isin(index)
        assert result.all()

        result = index.isin(list(index))
        assert result.all()

        assert_almost_equal(index.isin([index[2], 5]),
                            np.array([False, False, True, False]))

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(10, 10,
                                    data_gen_f=lambda *args, **kwargs: randn(),
                                    r_idx_type='i', c_idx_type='dt')
        cols = df.columns.join(df.index, how='outer')
        joined = cols.join(df.columns)
        assert cols.dtype == np.dtype('O')
        assert cols.dtype == joined.dtype
        tm.assert_numpy_array_equal(cols.values, joined.values)

    def test_join_self(self, join_type):
        index = date_range('1/1/2000', periods=10)
        joined = index.join(index, how=join_type)
        assert index is joined

    def assert_index_parameters(self, index):
        assert index.freq == '40960N'
        assert index.inferred_freq == '40960N'

    def test_ns_index(self):
        nsamples = 400
        ns = int(1e9 / 24414)
        dtstart = np.datetime64('2012-09-20T00:00:00')

        dt = dtstart + np.arange(nsamples) * np.timedelta64(ns, 'ns')
        freq = ns * offsets.Nano()
        index = pd.DatetimeIndex(dt, freq=freq, name='time')
        self.assert_index_parameters(index)

        new_index = pd.DatetimeIndex(start=index[0], end=index[-1],
                                     freq=index.freq)
        self.assert_index_parameters(new_index)

    def test_join_with_period_index(self, join_type):
        df = tm.makeCustomDataframe(
            10, 10, data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type='p', r_idx_type='dt')
        s = df.iloc[:5, 0]

        with tm.assert_raises_regex(ValueError,
                                    'can only call with other '
                                    'PeriodIndex-ed objects'):
            df.columns.join(s.index, how=join_type)

    def test_factorize(self):
        idx1 = DatetimeIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                              '2014-03', '2014-03'])

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])

        arr, idx = idx1.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # tz must be preserved
        idx1 = idx1.tz_localize('Asia/Tokyo')
        exp_idx = exp_idx.tz_localize('Asia/Tokyo')

        arr, idx = idx1.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        idx2 = pd.DatetimeIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                                 '2014-03', '2014-01'])

        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])
        arr, idx = idx2.factorize(sort=True)
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        exp_arr = np.array([0, 0, 1, 2, 0, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-03', '2014-02', '2014-01'])
        arr, idx = idx2.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # freq must be preserved
        idx3 = date_range('2000-01', periods=4, freq='M', tz='Asia/Tokyo')
        exp_arr = np.array([0, 1, 2, 3], dtype=np.intp)
        arr, idx = idx3.factorize()
        tm.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, idx3)

    def test_factorize_tz(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH#13750
        base = pd.date_range('2016-11-05', freq='H', periods=100, tz=tz)
        idx = base.repeat(5)

        exp_arr = np.arange(100, dtype=np.intp).repeat(5)

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, exp_arr)
            tm.assert_index_equal(res, base)

    def test_factorize_dst(self):
        # GH 13750
        idx = pd.date_range('2016-11-06', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

        idx = pd.date_range('2016-06-13', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            tm.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

    @pytest.mark.parametrize('arr, expected', [
        (pd.DatetimeIndex(['2017', '2017']), pd.DatetimeIndex(['2017'])),
        (pd.DatetimeIndex(['2017', '2017'], tz='US/Eastern'),
         pd.DatetimeIndex(['2017'], tz='US/Eastern')),
    ])
    def test_unique(self, arr, expected):
        result = arr.unique()
        tm.assert_index_equal(result, expected)
