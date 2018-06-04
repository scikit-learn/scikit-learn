# -*- coding: utf-8 -*-

from __future__ import print_function

from datetime import datetime, time

import pytest

from numpy import nan
from numpy.random import randn
import numpy as np

from pandas import (DataFrame, Series, Index,
                    Timestamp, DatetimeIndex, MultiIndex,
                    to_datetime, date_range, period_range)
import pandas as pd
import pandas.tseries.offsets as offsets

from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_index_equal,
                                 assert_raises_regex)

import pandas.util.testing as tm
from pandas.compat import product

from pandas.tests.frame.common import TestData


class TestDataFrameTimeSeriesMethods(TestData):

    def test_diff(self):
        the_diff = self.tsframe.diff(1)

        assert_series_equal(the_diff['A'],
                            self.tsframe['A'] - self.tsframe['A'].shift(1))

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = DataFrame({'s': s}).diff()
        assert rs.s[1] == 1

        # mixed numeric
        tf = self.tsframe.astype('float32')
        the_diff = tf.diff(1)
        assert_series_equal(the_diff['A'],
                            tf['A'] - tf['A'].shift(1))

        # issue 10907
        df = pd.DataFrame({'y': pd.Series([2]), 'z': pd.Series([3])})
        df.insert(0, 'x', 1)
        result = df.diff(axis=1)
        expected = pd.DataFrame({'x': np.nan, 'y': pd.Series(
            1), 'z': pd.Series(1)}).astype('float64')
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('tz', [None, 'UTC'])
    def test_diff_datetime_axis0(self, tz):
        # GH 18578
        df = DataFrame({0: date_range('2010', freq='D', periods=2, tz=tz),
                        1: date_range('2010', freq='D', periods=2, tz=tz)})

        result = df.diff(axis=0)
        expected = DataFrame({0: pd.TimedeltaIndex(['NaT', '1 days']),
                              1: pd.TimedeltaIndex(['NaT', '1 days'])})
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('tz', [None, 'UTC'])
    def test_diff_datetime_axis1(self, tz):
        # GH 18578
        df = DataFrame({0: date_range('2010', freq='D', periods=2, tz=tz),
                        1: date_range('2010', freq='D', periods=2, tz=tz)})
        if tz is None:
            result = df.diff(axis=1)
            expected = DataFrame({0: pd.TimedeltaIndex(['NaT', 'NaT']),
                                  1: pd.TimedeltaIndex(['0 days',
                                                        '0 days'])})
            assert_frame_equal(result, expected)
        else:
            with pytest.raises(NotImplementedError):
                result = df.diff(axis=1)

    def test_diff_timedelta(self):
        # GH 4533
        df = DataFrame(dict(time=[Timestamp('20130101 9:01'),
                                  Timestamp('20130101 9:02')],
                            value=[1.0, 2.0]))

        res = df.diff()
        exp = DataFrame([[pd.NaT, np.nan],
                         [pd.Timedelta('00:01:00'), 1]],
                        columns=['time', 'value'])
        assert_frame_equal(res, exp)

    def test_diff_mixed_dtype(self):
        df = DataFrame(np.random.randn(5, 3))
        df['A'] = np.array([1, 2, 3, 4, 5], dtype=object)

        result = df.diff()
        assert result[0].dtype == np.float64

    def test_diff_neg_n(self):
        rs = self.tsframe.diff(-1)
        xp = self.tsframe - self.tsframe.shift(-1)
        assert_frame_equal(rs, xp)

    def test_diff_float_n(self):
        rs = self.tsframe.diff(1.)
        xp = self.tsframe.diff(1)
        assert_frame_equal(rs, xp)

    def test_diff_axis(self):
        # GH 9727
        df = DataFrame([[1., 2.], [3., 4.]])
        assert_frame_equal(df.diff(axis=1), DataFrame(
            [[np.nan, 1.], [np.nan, 1.]]))
        assert_frame_equal(df.diff(axis=0), DataFrame(
            [[np.nan, np.nan], [2., 2.]]))

    def test_pct_change(self):
        rs = self.tsframe.pct_change(fill_method=None)
        assert_frame_equal(rs, self.tsframe / self.tsframe.shift(1) - 1)

        rs = self.tsframe.pct_change(2)
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs, filled / filled.shift(2) - 1)

        rs = self.tsframe.pct_change(fill_method='bfill', limit=1)
        filled = self.tsframe.fillna(method='bfill', limit=1)
        assert_frame_equal(rs, filled / filled.shift(1) - 1)

        rs = self.tsframe.pct_change(freq='5D')
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs,
                           (filled / filled.shift(freq='5D') - 1)
                           .reindex_like(filled))

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        df = DataFrame({'a': s, 'b': s})

        chg = df.pct_change()
        expected = Series([np.nan, 0.5, 0., 2.5 / 1.5 - 1, .2])
        edf = DataFrame({'a': expected, 'b': expected})
        assert_frame_equal(chg, edf)

    @pytest.mark.parametrize("freq, periods, fill_method, limit",
                             [('5B', 5, None, None),
                              ('3B', 3, None, None),
                              ('3B', 3, 'bfill', None),
                              ('7B', 7, 'pad', 1),
                              ('7B', 7, 'bfill', 3),
                              ('14B', 14, None, None)])
    def test_pct_change_periods_freq(self, freq, periods, fill_method, limit):
        # GH 7292
        rs_freq = self.tsframe.pct_change(freq=freq,
                                          fill_method=fill_method,
                                          limit=limit)
        rs_periods = self.tsframe.pct_change(periods,
                                             fill_method=fill_method,
                                             limit=limit)
        assert_frame_equal(rs_freq, rs_periods)

        empty_ts = DataFrame(index=self.tsframe.index,
                             columns=self.tsframe.columns)
        rs_freq = empty_ts.pct_change(freq=freq,
                                      fill_method=fill_method,
                                      limit=limit)
        rs_periods = empty_ts.pct_change(periods,
                                         fill_method=fill_method,
                                         limit=limit)
        assert_frame_equal(rs_freq, rs_periods)

    def test_frame_ctor_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50', freq='10s')
        dates = np.asarray(rng)

        df = DataFrame({'A': np.random.randn(len(rng)), 'B': dates})
        assert np.issubdtype(df['B'].dtype, np.dtype('M8[ns]'))

    def test_frame_add_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50', freq='10s')
        df = DataFrame(index=np.arange(len(rng)))

        df['A'] = rng
        assert np.issubdtype(df['A'].dtype, np.dtype('M8[ns]'))

    def test_frame_datetime64_pre1900_repr(self):
        df = DataFrame({'year': date_range('1/1/1700', periods=50,
                                           freq='A-DEC')})
        # it works!
        repr(df)

    def test_frame_add_datetime64_col_other_units(self):
        n = 100

        units = ['h', 'm', 's', 'ms', 'D', 'M', 'Y']

        ns_dtype = np.dtype('M8[ns]')

        for unit in units:
            dtype = np.dtype('M8[%s]' % unit)
            vals = np.arange(n, dtype=np.int64).view(dtype)

            df = DataFrame({'ints': np.arange(n)}, index=np.arange(n))
            df[unit] = vals

            ex_vals = to_datetime(vals.astype('O')).values

            assert df[unit].dtype == ns_dtype
            assert (df[unit].values == ex_vals).all()

        # Test insertion into existing datetime64 column
        df = DataFrame({'ints': np.arange(n)}, index=np.arange(n))
        df['dates'] = np.arange(n, dtype=np.int64).view(ns_dtype)

        for unit in units:
            dtype = np.dtype('M8[%s]' % unit)
            vals = np.arange(n, dtype=np.int64).view(dtype)

            tmp = df.copy()

            tmp['dates'] = vals
            ex_vals = to_datetime(vals.astype('O')).values

            assert (tmp['dates'].values == ex_vals).all()

    def test_shift(self):
        # naive shift
        shiftedFrame = self.tsframe.shift(5)
        tm.assert_index_equal(shiftedFrame.index, self.tsframe.index)

        shiftedSeries = self.tsframe['A'].shift(5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        shiftedFrame = self.tsframe.shift(-5)
        tm.assert_index_equal(shiftedFrame.index, self.tsframe.index)

        shiftedSeries = self.tsframe['A'].shift(-5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        # shift by 0
        unshifted = self.tsframe.shift(0)
        assert_frame_equal(unshifted, self.tsframe)

        # shift by DateOffset
        shiftedFrame = self.tsframe.shift(5, freq=offsets.BDay())
        assert len(shiftedFrame) == len(self.tsframe)

        shiftedFrame2 = self.tsframe.shift(5, freq='B')
        assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = self.tsframe.index[0]
        shifted_d = d + offsets.BDay(5)
        assert_series_equal(self.tsframe.xs(d),
                            shiftedFrame.xs(shifted_d), check_names=False)

        # shift int frame
        int_shifted = self.intframe.shift(1)  # noqa

        # Shifting with PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_index_equal(shifted.index, ps.index)
        tm.assert_index_equal(unshifted.index, ps.index)
        tm.assert_numpy_array_equal(unshifted.iloc[:, 0].dropna().values,
                                    ps.iloc[:-1, 0].values)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, offsets.BDay())
        assert_frame_equal(shifted2, shifted3)
        assert_frame_equal(ps, shifted2.shift(-1, 'B'))

        tm.assert_raises_regex(ValueError,
                               'does not match PeriodIndex freq',
                               ps.shift, freq='D')

        # shift other axis
        # GH 6371
        df = DataFrame(np.random.rand(10, 5))
        expected = pd.concat([DataFrame(np.nan, index=df.index,
                                        columns=[0]),
                              df.iloc[:, 0:-1]],
                             ignore_index=True, axis=1)
        result = df.shift(1, axis=1)
        assert_frame_equal(result, expected)

        # shift named axis
        df = DataFrame(np.random.rand(10, 5))
        expected = pd.concat([DataFrame(np.nan, index=df.index,
                                        columns=[0]),
                              df.iloc[:, 0:-1]],
                             ignore_index=True, axis=1)
        result = df.shift(1, axis='columns')
        assert_frame_equal(result, expected)

    def test_shift_bool(self):
        df = DataFrame({'high': [True, False],
                        'low': [False, False]})
        rs = df.shift(1)
        xp = DataFrame(np.array([[np.nan, np.nan],
                                 [True, False]], dtype=object),
                       columns=['high', 'low'])
        assert_frame_equal(rs, xp)

    def test_shift_categorical(self):
        # GH 9416
        s1 = pd.Series(['a', 'b', 'c'], dtype='category')
        s2 = pd.Series(['A', 'B', 'C'], dtype='category')
        df = DataFrame({'one': s1, 'two': s2})
        rs = df.shift(1)
        xp = DataFrame({'one': s1.shift(1), 'two': s2.shift(1)})
        assert_frame_equal(rs, xp)

    def test_shift_empty(self):
        # Regression test for #8019
        df = DataFrame({'foo': []})
        rs = df.shift(-1)

        assert_frame_equal(df, rs)

    def test_shift_duplicate_columns(self):
        # GH 9092; verify that position-based shifting works
        # in the presence of duplicate columns
        column_lists = [list(range(5)), [1] * 5, [1, 1, 2, 2, 1]]
        data = np.random.randn(20, 5)

        shifted = []
        for columns in column_lists:
            df = pd.DataFrame(data.copy(), columns=columns)
            for s in range(5):
                df.iloc[:, s] = df.iloc[:, s].shift(s + 1)
            df.columns = range(5)
            shifted.append(df)

        # sanity check the base case
        nulls = shifted[0].isna().sum()
        assert_series_equal(nulls, Series(range(1, 6), dtype='int64'))

        # check all answers are the same
        assert_frame_equal(shifted[0], shifted[1])
        assert_frame_equal(shifted[0], shifted[2])

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_frame_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=offsets.BDay())
        assert_frame_equal(shifted, shifted3)

        tm.assert_raises_regex(
            ValueError, 'does not match', ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.tsframe.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(self.tsframe, unshifted)

        shifted2 = self.tsframe.tshift(freq=self.tsframe.index.freq)
        assert_frame_equal(shifted, shifted2)

        inferred_ts = DataFrame(self.tsframe.values,
                                Index(np.asarray(self.tsframe.index)),
                                columns=self.tsframe.columns)
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_frame_equal(shifted, self.tsframe.tshift(1))
        assert_frame_equal(unshifted, inferred_ts)

        no_freq = self.tsframe.iloc[[0, 5, 7], :]
        pytest.raises(ValueError, no_freq.tshift)

    def test_truncate(self):
        ts = self.tsframe[::3]

        start, end = self.tsframe.index[3], self.tsframe.index[6]

        start_missing = self.tsframe.index[2]
        end_missing = self.tsframe.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_frame_equal(truncated, expected)

        pytest.raises(ValueError, ts.truncate,
                      before=ts.index[-1] - 1,
                      after=ts.index[0] + 1)

    def test_truncate_copy(self):
        index = self.tsframe.index
        truncated = self.tsframe.truncate(index[5], index[10])
        truncated.values[:] = 5.
        assert not (self.tsframe.values[5:11] == 5).any()

    def test_truncate_nonsortedindex(self):
        # GH 17935

        df = pd.DataFrame({'A': ['a', 'b', 'c', 'd', 'e']},
                          index=[5, 3, 2, 9, 0])
        with tm.assert_raises_regex(ValueError,
                                    'truncate requires a sorted index'):
            df.truncate(before=3, after=9)

        rng = pd.date_range('2011-01-01', '2012-01-01', freq='W')
        ts = pd.DataFrame({'A': np.random.randn(len(rng)),
                           'B': np.random.randn(len(rng))},
                          index=rng)
        with tm.assert_raises_regex(ValueError,
                                    'truncate requires a sorted index'):
            ts.sort_values('A', ascending=False).truncate(before='2011-11',
                                                          after='2011-12')

        df = pd.DataFrame({3: np.random.randn(5),
                           20: np.random.randn(5),
                           2: np.random.randn(5),
                           0: np.random.randn(5)},
                          columns=[3, 20, 2, 0])
        with tm.assert_raises_regex(ValueError,
                                    'truncate requires a sorted index'):
            df.truncate(before=2, after=20, axis=1)

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(offsets.BMonthEnd())
        rule_monthly = self.tsframe.asfreq('BM')

        tm.assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('B', method='pad')  # noqa
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq('B', method='pad')  # noqa

        # test does not blow up on length-0 DataFrame
        zero_length = self.tsframe.reindex([])
        result = zero_length.asfreq('BM')
        assert result is not zero_length

    def test_asfreq_datetimeindex(self):
        df = DataFrame({'A': [1, 2, 3]},
                       index=[datetime(2011, 11, 1), datetime(2011, 11, 2),
                              datetime(2011, 11, 3)])
        df = df.asfreq('B')
        assert isinstance(df.index, DatetimeIndex)

        ts = df['A'].asfreq('B')
        assert isinstance(ts.index, DatetimeIndex)

    def test_asfreq_fillvalue(self):
        # test for fill value during upsampling, related to issue 3715

        # setup
        rng = pd.date_range('1/1/2016', periods=10, freq='2S')
        ts = pd.Series(np.arange(len(rng)), index=rng)
        df = pd.DataFrame({'one': ts})

        # insert pre-existing missing value
        df.loc['2016-01-01 00:00:08', 'one'] = None

        actual_df = df.asfreq(freq='1S', fill_value=9.0)
        expected_df = df.asfreq(freq='1S').fillna(9.0)
        expected_df.loc['2016-01-01 00:00:08', 'one'] = None
        assert_frame_equal(expected_df, actual_df)

        expected_series = ts.asfreq(freq='1S').fillna(9.0)
        actual_series = ts.asfreq(freq='1S', fill_value=9.0)
        assert_series_equal(expected_series, actual_series)

    def test_first_last_valid(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan
        mat[-5:] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        index = frame.first_valid_index()

        assert index == frame.index[5]

        index = frame.last_valid_index()
        assert index == frame.index[-6]

        # GH12800
        empty = DataFrame()
        assert empty.last_valid_index() is None
        assert empty.first_valid_index() is None

        # GH17400: no valid entries
        frame[:] = nan
        assert frame.last_valid_index() is None
        assert frame.first_valid_index() is None

        # GH20499: its preserves freq with holes
        frame.index = date_range("20110101", periods=N, freq="B")
        frame.iloc[1] = 1
        frame.iloc[-2] = 1
        assert frame.first_valid_index() == frame.index[1]
        assert frame.last_valid_index() == frame.index[-2]
        assert frame.first_valid_index().freq == frame.index.freq
        assert frame.last_valid_index().freq == frame.index.freq

    def test_first_subset(self):
        ts = tm.makeTimeDataFrame(freq='12h')
        result = ts.first('10d')
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(freq='D')
        result = ts.first('10d')
        assert len(result) == 10

        result = ts.first('3M')
        expected = ts[:'3/31/2000']
        assert_frame_equal(result, expected)

        result = ts.first('21D')
        expected = ts[:21]
        assert_frame_equal(result, expected)

        result = ts[:0].first('3M')
        assert_frame_equal(result, ts[:0])

    def test_first_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.first('1D')

    def test_last_subset(self):
        ts = tm.makeTimeDataFrame(freq='12h')
        result = ts.last('10d')
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(nper=30, freq='D')
        result = ts.last('10d')
        assert len(result) == 10

        result = ts.last('21D')
        expected = ts['2000-01-10':]
        assert_frame_equal(result, expected)

        result = ts.last('21D')
        expected = ts[-21:]
        assert_frame_equal(result, expected)

        result = ts[:0].last('3M')
        assert_frame_equal(result, ts[:0])

    def test_last_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.last('1D')

    def test_at_time(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        rs = ts.at_time(rng[1])
        assert (rs.index.hour == rng[1].hour).all()
        assert (rs.index.minute == rng[1].minute).all()
        assert (rs.index.second == rng[1].second).all()

        result = ts.at_time('9:30')
        expected = ts.at_time(time(9, 30))
        assert_frame_equal(result, expected)

        result = ts.loc[time(9, 30)]
        expected = ts.loc[(rng.hour == 9) & (rng.minute == 30)]

        assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range('1/1/2000', '1/31/2000')
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts.at_time(time(0, 0))
        assert_frame_equal(result, ts)

        # time doesn't exist
        rng = date_range('1/1/2012', freq='23Min', periods=384)
        ts = DataFrame(np.random.randn(len(rng), 2), rng)
        rs = ts.at_time('16:00')
        assert len(rs) == 0

    def test_at_time_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.at_time('00:00')

    def test_between_time(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = 13 * 4 + 1
            if not inc_start:
                exp_len -= 5
            if not inc_end:
                exp_len -= 4

            assert len(filtered) == exp_len
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    assert t >= stime
                else:
                    assert t > stime

                if inc_end:
                    assert t <= etime
                else:
                    assert t < etime

        result = ts.between_time('00:00', '01:00')
        expected = ts.between_time(stime, etime)
        assert_frame_equal(result, expected)

        # across midnight
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = (12 * 11 + 1) * 4 + 1
            if not inc_start:
                exp_len -= 4
            if not inc_end:
                exp_len -= 4

            assert len(filtered) == exp_len
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    assert (t >= stime) or (t <= etime)
                else:
                    assert (t > stime) or (t <= etime)

                if inc_end:
                    assert (t <= etime) or (t >= stime)
                else:
                    assert (t < etime) or (t >= stime)

    def test_between_time_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.between_time(start_time='00:00', end_time='12:00')

    def test_operation_on_NaT(self):
        # Both NaT and Timestamp are in DataFrame.
        df = pd.DataFrame({'foo': [pd.NaT, pd.NaT,
                                   pd.Timestamp('2012-05-01')]})

        res = df.min()
        exp = pd.Series([pd.Timestamp('2012-05-01')], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.Timestamp('2012-05-01')], index=["foo"])
        tm.assert_series_equal(res, exp)

        # GH12941, only NaTs are in DataFrame.
        df = pd.DataFrame({'foo': [pd.NaT, pd.NaT]})

        res = df.min()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

    def test_datetime_assignment_with_NaT_and_diff_time_units(self):
        # GH 7492
        data_ns = np.array([1, 'nat'], dtype='datetime64[ns]')
        result = pd.Series(data_ns).to_frame()
        result['new'] = data_ns
        expected = pd.DataFrame({0: [1, None],
                                 'new': [1, None]}, dtype='datetime64[ns]')
        tm.assert_frame_equal(result, expected)
        # OutOfBoundsDatetime error shouldn't occur
        data_s = np.array([1, 'nat'], dtype='datetime64[s]')
        result['new'] = data_s
        expected = pd.DataFrame({0: [1, None],
                                 'new': [1e9, None]}, dtype='datetime64[ns]')
        tm.assert_frame_equal(result, expected)

    def test_frame_to_period(self):
        K = 5
        from pandas.core.indexes.period import period_range

        dr = date_range('1/1/2000', '1/1/2001')
        pr = period_range('1/1/2000', '1/1/2001')
        df = DataFrame(randn(len(dr), K), index=dr)
        df['mix'] = 'a'

        pts = df.to_period()
        exp = df.copy()
        exp.index = pr
        assert_frame_equal(pts, exp)

        pts = df.to_period('M')
        tm.assert_index_equal(pts.index, exp.index.asfreq('M'))

        df = df.T
        pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = pr
        assert_frame_equal(pts, exp)

        pts = df.to_period('M', axis=1)
        tm.assert_index_equal(pts.columns, exp.columns.asfreq('M'))

        pytest.raises(ValueError, df.to_period, axis=2)

    @pytest.mark.parametrize("fn", ['tz_localize', 'tz_convert'])
    def test_tz_convert_and_localize(self, fn):
        l0 = date_range('20140701', periods=5, freq='D')

        # TODO: l1 should be a PeriodIndex for testing
        #       after GH2106 is addressed
        with pytest.raises(NotImplementedError):
            period_range('20140701', periods=1).tz_convert('UTC')
        with pytest.raises(NotImplementedError):
            period_range('20140701', periods=1).tz_localize('UTC')
        # l1 = period_range('20140701', periods=5, freq='D')
        l1 = date_range('20140701', periods=5, freq='D')

        int_idx = Index(range(5))

        if fn == 'tz_convert':
            l0 = l0.tz_localize('UTC')
            l1 = l1.tz_localize('UTC')

        for idx in [l0, l1]:

            l0_expected = getattr(idx, fn)('US/Pacific')
            l1_expected = getattr(idx, fn)('US/Pacific')

            df1 = DataFrame(np.ones(5), index=l0)
            df1 = getattr(df1, fn)('US/Pacific')
            assert_index_equal(df1.index, l0_expected)

            # MultiIndex
            # GH7846
            df2 = DataFrame(np.ones(5), MultiIndex.from_arrays([l0, l1]))

            df3 = getattr(df2, fn)('US/Pacific', level=0)
            assert not df3.index.levels[0].equals(l0)
            assert_index_equal(df3.index.levels[0], l0_expected)
            assert_index_equal(df3.index.levels[1], l1)
            assert not df3.index.levels[1].equals(l1_expected)

            df3 = getattr(df2, fn)('US/Pacific', level=1)
            assert_index_equal(df3.index.levels[0], l0)
            assert not df3.index.levels[0].equals(l0_expected)
            assert_index_equal(df3.index.levels[1], l1_expected)
            assert not df3.index.levels[1].equals(l1)

            df4 = DataFrame(np.ones(5),
                            MultiIndex.from_arrays([int_idx, l0]))

            # TODO: untested
            df5 = getattr(df4, fn)('US/Pacific', level=1)  # noqa

            assert_index_equal(df3.index.levels[0], l0)
            assert not df3.index.levels[0].equals(l0_expected)
            assert_index_equal(df3.index.levels[1], l1_expected)
            assert not df3.index.levels[1].equals(l1)

        # Bad Inputs

        # Not DatetimeIndex / PeriodIndex
        with assert_raises_regex(TypeError, 'DatetimeIndex'):
            df = DataFrame(index=int_idx)
            df = getattr(df, fn)('US/Pacific')

        # Not DatetimeIndex / PeriodIndex
        with assert_raises_regex(TypeError, 'DatetimeIndex'):
            df = DataFrame(np.ones(5),
                           MultiIndex.from_arrays([int_idx, l0]))
            df = getattr(df, fn)('US/Pacific', level=0)

        # Invalid level
        with assert_raises_regex(ValueError, 'not valid'):
            df = DataFrame(index=l0)
            df = getattr(df, fn)('US/Pacific', level=1)
