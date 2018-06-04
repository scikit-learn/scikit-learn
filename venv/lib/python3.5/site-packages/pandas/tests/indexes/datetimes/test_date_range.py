"""
test date_range, bdate_range construction from the convenience range functions
"""

import pytest

import numpy as np
import pytz
from pytz import timezone
from datetime import datetime, timedelta, time

import pandas as pd
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas import compat
from pandas import date_range, bdate_range, offsets, DatetimeIndex, Timestamp
from pandas.tseries.offsets import (generate_range, CDay, BDay, DateOffset,
                                    MonthEnd, prefix_mapping)

from pandas.tests.series.common import TestData

START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)


class TestTimestampEquivDateRange(object):
    # Older tests in TestTimeSeries constructed their `stamp` objects
    # using `date_range` instead of the `Timestamp` constructor.
    # TestTimestampEquivDateRange checks that these are equivalent in the
    # pertinent cases.

    def test_date_range_timestamp_equiv(self):
        rng = date_range('20090415', '20090519', tz='US/Eastern')
        stamp = rng[0]

        ts = Timestamp('20090415', tz='US/Eastern', freq='D')
        assert ts == stamp

    def test_date_range_timestamp_equiv_dateutil(self):
        rng = date_range('20090415', '20090519', tz='dateutil/US/Eastern')
        stamp = rng[0]

        ts = Timestamp('20090415', tz='dateutil/US/Eastern', freq='D')
        assert ts == stamp

    def test_date_range_timestamp_equiv_explicit_pytz(self):
        rng = date_range('20090415', '20090519',
                         tz=pytz.timezone('US/Eastern'))
        stamp = rng[0]

        ts = Timestamp('20090415', tz=pytz.timezone('US/Eastern'), freq='D')
        assert ts == stamp

    @td.skip_if_windows_python_3
    def test_date_range_timestamp_equiv_explicit_dateutil(self):
        from pandas._libs.tslibs.timezones import dateutil_gettz as gettz

        rng = date_range('20090415', '20090519', tz=gettz('US/Eastern'))
        stamp = rng[0]

        ts = Timestamp('20090415', tz=gettz('US/Eastern'), freq='D')
        assert ts == stamp

    def test_date_range_timestamp_equiv_from_datetime_instance(self):
        datetime_instance = datetime(2014, 3, 4)
        # build a timestamp with a frequency, since then it supports
        # addition/subtraction of integers
        timestamp_instance = date_range(datetime_instance, periods=1,
                                        freq='D')[0]

        ts = Timestamp(datetime_instance, freq='D')
        assert ts == timestamp_instance

    def test_date_range_timestamp_equiv_preserve_frequency(self):
        timestamp_instance = date_range('2014-03-05', periods=1, freq='D')[0]
        ts = Timestamp('2014-03-05', freq='D')

        assert timestamp_instance == ts


class TestDateRanges(TestData):

    def test_date_range_gen_error(self):
        rng = date_range('1/1/2000 00:00', '1/1/2000 00:18', freq='5min')
        assert len(rng) == 4

    @pytest.mark.parametrize("freq", ["AS", "YS"])
    def test_begin_year_alias(self, freq):
        # see gh-9313
        rng = date_range("1/1/2013", "7/1/2017", freq=freq)
        exp = pd.DatetimeIndex(["2013-01-01", "2014-01-01",
                                "2015-01-01", "2016-01-01",
                                "2017-01-01"], freq=freq)
        tm.assert_index_equal(rng, exp)

    @pytest.mark.parametrize("freq", ["A", "Y"])
    def test_end_year_alias(self, freq):
        # see gh-9313
        rng = date_range("1/1/2013", "7/1/2017", freq=freq)
        exp = pd.DatetimeIndex(["2013-12-31", "2014-12-31",
                                "2015-12-31", "2016-12-31"], freq=freq)
        tm.assert_index_equal(rng, exp)

    @pytest.mark.parametrize("freq", ["BA", "BY"])
    def test_business_end_year_alias(self, freq):
        # see gh-9313
        rng = date_range("1/1/2013", "7/1/2017", freq=freq)
        exp = pd.DatetimeIndex(["2013-12-31", "2014-12-31",
                                "2015-12-31", "2016-12-30"], freq=freq)
        tm.assert_index_equal(rng, exp)

    def test_date_range_negative_freq(self):
        # GH 11018
        rng = date_range('2011-12-31', freq='-2A', periods=3)
        exp = pd.DatetimeIndex(['2011-12-31', '2009-12-31',
                                '2007-12-31'], freq='-2A')
        tm.assert_index_equal(rng, exp)
        assert rng.freq == '-2A'

        rng = date_range('2011-01-31', freq='-2M', periods=3)
        exp = pd.DatetimeIndex(['2011-01-31', '2010-11-30',
                                '2010-09-30'], freq='-2M')
        tm.assert_index_equal(rng, exp)
        assert rng.freq == '-2M'

    def test_date_range_bms_bug(self):
        # #1645
        rng = date_range('1/1/2000', periods=10, freq='BMS')

        ex_first = Timestamp('2000-01-03')
        assert rng[0] == ex_first

    def test_date_range_normalize(self):
        snap = datetime.today()
        n = 50

        rng = date_range(snap, periods=n, normalize=False, freq='2D')

        offset = timedelta(2)
        values = DatetimeIndex([snap + i * offset for i in range(n)])

        tm.assert_index_equal(rng, values)

        rng = date_range('1/1/2000 08:15', periods=n, normalize=False,
                         freq='B')
        the_time = time(8, 15)
        for val in rng:
            assert val.time() == the_time

    def test_date_range_fy5252(self):
        dr = date_range(start="2013-01-01", periods=2, freq=offsets.FY5253(
            startingMonth=1, weekday=3, variation="nearest"))
        assert dr[0] == Timestamp('2013-01-31')
        assert dr[1] == Timestamp('2014-01-30')

    def test_date_range_ambiguous_arguments(self):
        # #2538
        start = datetime(2011, 1, 1, 5, 3, 40)
        end = datetime(2011, 1, 1, 8, 9, 40)

        msg = ('Of the four parameters: start, end, periods, and '
               'freq, exactly three must be specified')
        with tm.assert_raises_regex(ValueError, msg):
            date_range(start, end, periods=10, freq='s')

    def test_date_range_convenience_periods(self):
        # GH 20808
        result = date_range('2018-04-24', '2018-04-27', periods=3)
        expected = DatetimeIndex(['2018-04-24 00:00:00',
                                  '2018-04-25 12:00:00',
                                  '2018-04-27 00:00:00'], freq=None)

        tm.assert_index_equal(result, expected)

        # Test if spacing remains linear if tz changes to dst in range
        result = date_range('2018-04-01 01:00:00',
                            '2018-04-01 04:00:00',
                            tz='Australia/Sydney',
                            periods=3)
        expected = DatetimeIndex([Timestamp('2018-04-01 01:00:00+1100',
                                            tz='Australia/Sydney'),
                                  Timestamp('2018-04-01 02:00:00+1000',
                                            tz='Australia/Sydney'),
                                  Timestamp('2018-04-01 04:00:00+1000',
                                            tz='Australia/Sydney')])
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('start,end,result_tz', [
        ['20180101', '20180103', 'US/Eastern'],
        [datetime(2018, 1, 1), datetime(2018, 1, 3), 'US/Eastern'],
        [Timestamp('20180101'), Timestamp('20180103'), 'US/Eastern'],
        [Timestamp('20180101', tz='US/Eastern'),
         Timestamp('20180103', tz='US/Eastern'), 'US/Eastern'],
        [Timestamp('20180101', tz='US/Eastern'),
         Timestamp('20180103', tz='US/Eastern'), None]])
    def test_date_range_linspacing_tz(self, start, end, result_tz):
        # GH 20983
        result = date_range(start, end, periods=3, tz=result_tz)
        expected = date_range('20180101', periods=3, freq='D', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

    def test_date_range_businesshour(self):
        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00',
                             '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00',
                             '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00'],
                            freq='BH')
        rng = date_range('2014-07-04 09:00', '2014-07-04 16:00', freq='BH')
        tm.assert_index_equal(idx, rng)

        idx = DatetimeIndex(
            ['2014-07-04 16:00', '2014-07-07 09:00'], freq='BH')
        rng = date_range('2014-07-04 16:00', '2014-07-07 09:00', freq='BH')
        tm.assert_index_equal(idx, rng)

        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00',
                             '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00',
                             '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00',
                             '2014-07-07 09:00', '2014-07-07 10:00',
                             '2014-07-07 11:00',
                             '2014-07-07 12:00', '2014-07-07 13:00',
                             '2014-07-07 14:00',
                             '2014-07-07 15:00', '2014-07-07 16:00',
                             '2014-07-08 09:00', '2014-07-08 10:00',
                             '2014-07-08 11:00',
                             '2014-07-08 12:00', '2014-07-08 13:00',
                             '2014-07-08 14:00',
                             '2014-07-08 15:00', '2014-07-08 16:00'],
                            freq='BH')
        rng = date_range('2014-07-04 09:00', '2014-07-08 16:00', freq='BH')
        tm.assert_index_equal(idx, rng)

    def test_range_misspecified(self):
        # GH #1095
        msg = ('Of the four parameters: start, end, periods, and '
               'freq, exactly three must be specified')

        with tm.assert_raises_regex(ValueError, msg):
            date_range(start='1/1/2000')

        with tm.assert_raises_regex(ValueError, msg):
            date_range(end='1/1/2000')

        with tm.assert_raises_regex(ValueError, msg):
            date_range(periods=10)

        with tm.assert_raises_regex(ValueError, msg):
            date_range(start='1/1/2000', freq='H')

        with tm.assert_raises_regex(ValueError, msg):
            date_range(end='1/1/2000', freq='H')

        with tm.assert_raises_regex(ValueError, msg):
            date_range(periods=10, freq='H')

        with tm.assert_raises_regex(ValueError, msg):
            date_range()

    @pytest.mark.parametrize('f', [compat.long, int])
    def test_compat_replace(self, f):
        # https://github.com/statsmodels/statsmodels/issues/3349
        # replace should take ints/longs for compat
        result = date_range(Timestamp('1960-04-01 00:00:00', freq='QS-JAN'),
                            periods=f(76), freq='QS-JAN')
        assert len(result) == 76

    def test_catch_infinite_loop(self):
        offset = offsets.DateOffset(minute=5)
        # blow up, don't loop forever
        pytest.raises(Exception, date_range, datetime(2011, 11, 11),
                      datetime(2011, 11, 12), freq=offset)

    @pytest.mark.parametrize('periods', (1, 2))
    def test_wom_len(self, periods):
        # https://github.com/pandas-dev/pandas/issues/20517
        res = date_range(start='20110101', periods=periods, freq='WOM-1MON')
        assert len(res) == periods


class TestGenRangeGeneration(object):

    def test_generate(self):
        rng1 = list(generate_range(START, END, offset=BDay()))
        rng2 = list(generate_range(START, END, time_rule='B'))
        assert rng1 == rng2

    def test_generate_cday(self):
        rng1 = list(generate_range(START, END, offset=CDay()))
        rng2 = list(generate_range(START, END, time_rule='C'))
        assert rng1 == rng2

    def test_1(self):
        rng = list(generate_range(start=datetime(2009, 3, 25), periods=2))
        expected = [datetime(2009, 3, 25), datetime(2009, 3, 26)]
        assert rng == expected

    def test_2(self):
        rng = list(generate_range(start=datetime(2008, 1, 1),
                                  end=datetime(2008, 1, 3)))
        expected = [datetime(2008, 1, 1),
                    datetime(2008, 1, 2),
                    datetime(2008, 1, 3)]
        assert rng == expected

    def test_3(self):
        rng = list(generate_range(start=datetime(2008, 1, 5),
                                  end=datetime(2008, 1, 6)))
        expected = []
        assert rng == expected

    def test_precision_finer_than_offset(self):
        # GH 9907
        result1 = DatetimeIndex(start='2015-04-15 00:00:03',
                                end='2016-04-22 00:00:00', freq='Q')
        result2 = DatetimeIndex(start='2015-04-15 00:00:03',
                                end='2015-06-22 00:00:04', freq='W')
        expected1_list = ['2015-06-30 00:00:03', '2015-09-30 00:00:03',
                          '2015-12-31 00:00:03', '2016-03-31 00:00:03']
        expected2_list = ['2015-04-19 00:00:03', '2015-04-26 00:00:03',
                          '2015-05-03 00:00:03', '2015-05-10 00:00:03',
                          '2015-05-17 00:00:03', '2015-05-24 00:00:03',
                          '2015-05-31 00:00:03', '2015-06-07 00:00:03',
                          '2015-06-14 00:00:03', '2015-06-21 00:00:03']
        expected1 = DatetimeIndex(expected1_list, dtype='datetime64[ns]',
                                  freq='Q-DEC', tz=None)
        expected2 = DatetimeIndex(expected2_list, dtype='datetime64[ns]',
                                  freq='W-SUN', tz=None)
        tm.assert_index_equal(result1, expected1)
        tm.assert_index_equal(result2, expected2)

    dt1, dt2 = '2017-01-01', '2017-01-01'
    tz1, tz2 = 'US/Eastern', 'Europe/London'

    @pytest.mark.parametrize("start,end", [
        (pd.Timestamp(dt1, tz=tz1), pd.Timestamp(dt2)),
        (pd.Timestamp(dt1), pd.Timestamp(dt2, tz=tz2)),
        (pd.Timestamp(dt1, tz=tz1), pd.Timestamp(dt2, tz=tz2)),
        (pd.Timestamp(dt1, tz=tz2), pd.Timestamp(dt2, tz=tz1))
    ])
    def test_mismatching_tz_raises_err(self, start, end):
        # issue 18488
        with pytest.raises(TypeError):
            pd.date_range(start, end)
        with pytest.raises(TypeError):
            pd.DatetimeIndex(start, end, freq=BDay())


class TestBusinessDateRange(object):

    def test_constructor(self):
        bdate_range(START, END, freq=BDay())
        bdate_range(START, periods=20, freq=BDay())
        bdate_range(end=START, periods=20, freq=BDay())

        msg = 'periods must be a number, got B'
        with tm.assert_raises_regex(TypeError, msg):
            date_range('2011-1-1', '2012-1-1', 'B')

        with tm.assert_raises_regex(TypeError, msg):
            bdate_range('2011-1-1', '2012-1-1', 'B')

        msg = 'freq must be specified for bdate_range; use date_range instead'
        with tm.assert_raises_regex(TypeError, msg):
            bdate_range(START, END, periods=10, freq=None)

    def test_naive_aware_conflicts(self):
        naive = bdate_range(START, END, freq=BDay(), tz=None)
        aware = bdate_range(START, END, freq=BDay(), tz="Asia/Hong_Kong")

        msg = 'tz-naive.*tz-aware'
        with tm.assert_raises_regex(TypeError, msg):
            naive.join(aware)

        with tm.assert_raises_regex(TypeError, msg):
            aware.join(naive)

    def test_cached_range(self):
        DatetimeIndex._cached_range(START, END, freq=BDay())
        DatetimeIndex._cached_range(START, periods=20, freq=BDay())
        DatetimeIndex._cached_range(end=START, periods=20, freq=BDay())

        with tm.assert_raises_regex(TypeError, "freq"):
            DatetimeIndex._cached_range(START, END)

        with tm.assert_raises_regex(TypeError, "specify period"):
            DatetimeIndex._cached_range(START, freq=BDay())

        with tm.assert_raises_regex(TypeError, "specify period"):
            DatetimeIndex._cached_range(end=END, freq=BDay())

        with tm.assert_raises_regex(TypeError, "start or end"):
            DatetimeIndex._cached_range(periods=20, freq=BDay())

    def test_cached_range_bug(self):
        rng = date_range('2010-09-01 05:00:00', periods=50,
                         freq=DateOffset(hours=6))
        assert len(rng) == 50
        assert rng[0] == datetime(2010, 9, 1, 5)

    def test_timezone_comparaison_bug(self):
        # smoke test
        start = Timestamp('20130220 10:00', tz='US/Eastern')
        result = date_range(start, periods=2, tz='US/Eastern')
        assert len(result) == 2

    def test_timezone_comparaison_assert(self):
        start = Timestamp('20130220 10:00', tz='US/Eastern')
        msg = 'Inferred time zone not equal to passed time zone'
        with tm.assert_raises_regex(AssertionError, msg):
            date_range(start, periods=2, tz='Europe/Berlin')

    def test_misc(self):
        end = datetime(2009, 5, 13)
        dr = bdate_range(end=end, periods=20)
        firstDate = end - 19 * BDay()

        assert len(dr) == 20
        assert dr[0] == firstDate
        assert dr[-1] == end

    def test_date_parse_failure(self):
        badly_formed_date = '2007/100/1'

        with pytest.raises(ValueError):
            Timestamp(badly_formed_date)

        with pytest.raises(ValueError):
            bdate_range(start=badly_formed_date, periods=10)

        with pytest.raises(ValueError):
            bdate_range(end=badly_formed_date, periods=10)

        with pytest.raises(ValueError):
            bdate_range(badly_formed_date, badly_formed_date)

    def test_daterange_bug_456(self):
        # GH #456
        rng1 = bdate_range('12/5/2011', '12/5/2011')
        rng2 = bdate_range('12/2/2011', '12/5/2011')
        rng2.freq = BDay()

        result = rng1.union(rng2)
        assert isinstance(result, DatetimeIndex)

    def test_error_with_zero_monthends(self):
        msg = r'Offset <0 \* MonthEnds> did not increment date'
        with tm.assert_raises_regex(ValueError, msg):
            date_range('1/1/2000', '1/1/2001', freq=MonthEnd(0))

    def test_range_bug(self):
        # GH #770
        offset = DateOffset(months=3)
        result = date_range("2011-1-1", "2012-1-31", freq=offset)

        start = datetime(2011, 1, 1)
        expected = DatetimeIndex([start + i * offset for i in range(5)])
        tm.assert_index_equal(result, expected)

    def test_range_tz_pytz(self):
        # see gh-2906
        tz = timezone('US/Eastern')
        start = tz.localize(datetime(2011, 1, 1))
        end = tz.localize(datetime(2011, 1, 3))

        dr = date_range(start=start, periods=3)
        assert dr.tz.zone == tz.zone
        assert dr[0] == start
        assert dr[2] == end

        dr = date_range(end=end, periods=3)
        assert dr.tz.zone == tz.zone
        assert dr[0] == start
        assert dr[2] == end

        dr = date_range(start=start, end=end)
        assert dr.tz.zone == tz.zone
        assert dr[0] == start
        assert dr[2] == end

    def test_range_tz_dst_straddle_pytz(self):
        tz = timezone('US/Eastern')
        dates = [(tz.localize(datetime(2014, 3, 6)),
                  tz.localize(datetime(2014, 3, 12))),
                 (tz.localize(datetime(2013, 11, 1)),
                  tz.localize(datetime(2013, 11, 6)))]
        for (start, end) in dates:
            dr = date_range(start, end, freq='D')
            assert dr[0] == start
            assert dr[-1] == end
            assert np.all(dr.hour == 0)

            dr = date_range(start, end, freq='D', tz='US/Eastern')
            assert dr[0] == start
            assert dr[-1] == end
            assert np.all(dr.hour == 0)

            dr = date_range(start.replace(tzinfo=None), end.replace(
                tzinfo=None), freq='D', tz='US/Eastern')
            assert dr[0] == start
            assert dr[-1] == end
            assert np.all(dr.hour == 0)

    def test_range_tz_dateutil(self):
        # see gh-2906

        # Use maybe_get_tz to fix filename in tz under dateutil.
        from pandas._libs.tslibs.timezones import maybe_get_tz
        tz = lambda x: maybe_get_tz('dateutil/' + x)

        start = datetime(2011, 1, 1, tzinfo=tz('US/Eastern'))
        end = datetime(2011, 1, 3, tzinfo=tz('US/Eastern'))

        dr = date_range(start=start, periods=3)
        assert dr.tz == tz('US/Eastern')
        assert dr[0] == start
        assert dr[2] == end

        dr = date_range(end=end, periods=3)
        assert dr.tz == tz('US/Eastern')
        assert dr[0] == start
        assert dr[2] == end

        dr = date_range(start=start, end=end)
        assert dr.tz == tz('US/Eastern')
        assert dr[0] == start
        assert dr[2] == end

    @pytest.mark.parametrize('freq', ["1D", "3D", "2M", "7W", "3H", "A"])
    def test_range_closed(self, freq):
        begin = datetime(2011, 1, 1)
        end = datetime(2014, 1, 1)

        closed = date_range(begin, end, closed=None, freq=freq)
        left = date_range(begin, end, closed="left", freq=freq)
        right = date_range(begin, end, closed="right", freq=freq)
        expected_left = left
        expected_right = right

        if end == closed[-1]:
            expected_left = closed[:-1]
        if begin == closed[0]:
            expected_right = closed[1:]

        tm.assert_index_equal(expected_left, left)
        tm.assert_index_equal(expected_right, right)

    def test_range_closed_with_tz_aware_start_end(self):
        # GH12409, GH12684
        begin = Timestamp('2011/1/1', tz='US/Eastern')
        end = Timestamp('2014/1/1', tz='US/Eastern')

        for freq in ["1D", "3D", "2M", "7W", "3H", "A"]:
            closed = date_range(begin, end, closed=None, freq=freq)
            left = date_range(begin, end, closed="left", freq=freq)
            right = date_range(begin, end, closed="right", freq=freq)
            expected_left = left
            expected_right = right

            if end == closed[-1]:
                expected_left = closed[:-1]
            if begin == closed[0]:
                expected_right = closed[1:]

            tm.assert_index_equal(expected_left, left)
            tm.assert_index_equal(expected_right, right)

        begin = Timestamp('2011/1/1')
        end = Timestamp('2014/1/1')
        begintz = Timestamp('2011/1/1', tz='US/Eastern')
        endtz = Timestamp('2014/1/1', tz='US/Eastern')

        for freq in ["1D", "3D", "2M", "7W", "3H", "A"]:
            closed = date_range(begin, end, closed=None, freq=freq,
                                tz='US/Eastern')
            left = date_range(begin, end, closed="left", freq=freq,
                              tz='US/Eastern')
            right = date_range(begin, end, closed="right", freq=freq,
                               tz='US/Eastern')
            expected_left = left
            expected_right = right

            if endtz == closed[-1]:
                expected_left = closed[:-1]
            if begintz == closed[0]:
                expected_right = closed[1:]

            tm.assert_index_equal(expected_left, left)
            tm.assert_index_equal(expected_right, right)

    @pytest.mark.parametrize('closed', ['right', 'left', None])
    def test_range_closed_boundary(self, closed):
        # GH#11804
        right_boundary = date_range('2015-09-12', '2015-12-01',
                                    freq='QS-MAR', closed=closed)
        left_boundary = date_range('2015-09-01', '2015-09-12',
                                   freq='QS-MAR', closed=closed)
        both_boundary = date_range('2015-09-01', '2015-12-01',
                                   freq='QS-MAR', closed=closed)
        expected_right = expected_left = expected_both = both_boundary

        if closed == 'right':
            expected_left = both_boundary[1:]
        if closed == 'left':
            expected_right = both_boundary[:-1]
        if closed is None:
            expected_right = both_boundary[1:]
            expected_left = both_boundary[:-1]

        tm.assert_index_equal(right_boundary, expected_right)
        tm.assert_index_equal(left_boundary, expected_left)
        tm.assert_index_equal(both_boundary, expected_both)

    def test_years_only(self):
        # GH 6961
        dr = date_range('2014', '2015', freq='M')
        assert dr[0] == datetime(2014, 1, 31)
        assert dr[-1] == datetime(2014, 12, 31)

    def test_freq_divides_end_in_nanos(self):
        # GH 10885
        result_1 = date_range('2005-01-12 10:00', '2005-01-12 16:00',
                              freq='345min')
        result_2 = date_range('2005-01-13 10:00', '2005-01-13 16:00',
                              freq='345min')
        expected_1 = DatetimeIndex(['2005-01-12 10:00:00',
                                    '2005-01-12 15:45:00'],
                                   dtype='datetime64[ns]', freq='345T',
                                   tz=None)
        expected_2 = DatetimeIndex(['2005-01-13 10:00:00',
                                    '2005-01-13 15:45:00'],
                                   dtype='datetime64[ns]', freq='345T',
                                   tz=None)
        tm.assert_index_equal(result_1, expected_1)
        tm.assert_index_equal(result_2, expected_2)


class TestCustomDateRange(object):

    def test_constructor(self):
        bdate_range(START, END, freq=CDay())
        bdate_range(START, periods=20, freq=CDay())
        bdate_range(end=START, periods=20, freq=CDay())

        msg = 'periods must be a number, got C'
        with tm.assert_raises_regex(TypeError, msg):
            date_range('2011-1-1', '2012-1-1', 'C')

        with tm.assert_raises_regex(TypeError, msg):
            bdate_range('2011-1-1', '2012-1-1', 'C')

    def test_cached_range(self):
        DatetimeIndex._cached_range(START, END, freq=CDay())
        DatetimeIndex._cached_range(START, periods=20,
                                    freq=CDay())
        DatetimeIndex._cached_range(end=START, periods=20,
                                    freq=CDay())

        # with pytest.raises(TypeError):
        with tm.assert_raises_regex(TypeError, "freq"):
            DatetimeIndex._cached_range(START, END)

        # with pytest.raises(TypeError):
        with tm.assert_raises_regex(TypeError, "specify period"):
            DatetimeIndex._cached_range(START, freq=CDay())

        # with pytest.raises(TypeError):
        with tm.assert_raises_regex(TypeError, "specify period"):
            DatetimeIndex._cached_range(end=END, freq=CDay())

        # with pytest.raises(TypeError):
        with tm.assert_raises_regex(TypeError, "start or end"):
            DatetimeIndex._cached_range(periods=20, freq=CDay())

    def test_misc(self):
        end = datetime(2009, 5, 13)
        dr = bdate_range(end=end, periods=20, freq='C')
        firstDate = end - 19 * CDay()

        assert len(dr) == 20
        assert dr[0] == firstDate
        assert dr[-1] == end

    def test_daterange_bug_456(self):
        # GH #456
        rng1 = bdate_range('12/5/2011', '12/5/2011', freq='C')
        rng2 = bdate_range('12/2/2011', '12/5/2011', freq='C')
        rng2.freq = CDay()

        result = rng1.union(rng2)
        assert isinstance(result, DatetimeIndex)

    def test_cdaterange(self):
        result = bdate_range('2013-05-01', periods=3, freq='C')
        expected = DatetimeIndex(['2013-05-01', '2013-05-02', '2013-05-03'])
        tm.assert_index_equal(result, expected)

    def test_cdaterange_weekmask(self):
        result = bdate_range('2013-05-01', periods=3, freq='C',
                             weekmask='Sun Mon Tue Wed Thu')
        expected = DatetimeIndex(['2013-05-01', '2013-05-02', '2013-05-05'])
        tm.assert_index_equal(result, expected)

        # raise with non-custom freq
        msg = ('a custom frequency string is required when holidays or '
               'weekmask are passed, got frequency B')
        with tm.assert_raises_regex(ValueError, msg):
            bdate_range('2013-05-01', periods=3,
                        weekmask='Sun Mon Tue Wed Thu')

    def test_cdaterange_holidays(self):
        result = bdate_range('2013-05-01', periods=3, freq='C',
                             holidays=['2013-05-01'])
        expected = DatetimeIndex(['2013-05-02', '2013-05-03', '2013-05-06'])
        tm.assert_index_equal(result, expected)

        # raise with non-custom freq
        msg = ('a custom frequency string is required when holidays or '
               'weekmask are passed, got frequency B')
        with tm.assert_raises_regex(ValueError, msg):
            bdate_range('2013-05-01', periods=3, holidays=['2013-05-01'])

    def test_cdaterange_weekmask_and_holidays(self):
        result = bdate_range('2013-05-01', periods=3, freq='C',
                             weekmask='Sun Mon Tue Wed Thu',
                             holidays=['2013-05-01'])
        expected = DatetimeIndex(['2013-05-02', '2013-05-05', '2013-05-06'])
        tm.assert_index_equal(result, expected)

        # raise with non-custom freq
        msg = ('a custom frequency string is required when holidays or '
               'weekmask are passed, got frequency B')
        with tm.assert_raises_regex(ValueError, msg):
            bdate_range('2013-05-01', periods=3,
                        weekmask='Sun Mon Tue Wed Thu',
                        holidays=['2013-05-01'])

    @pytest.mark.parametrize('freq', [freq for freq in prefix_mapping
                                      if freq.startswith('C')])
    def test_all_custom_freq(self, freq):
        # should not raise
        bdate_range(START, END, freq=freq, weekmask='Mon Wed Fri',
                    holidays=['2009-03-14'])

        bad_freq = freq + 'FOO'
        msg = 'invalid custom frequency string: {freq}'
        with tm.assert_raises_regex(ValueError, msg.format(freq=bad_freq)):
            bdate_range(START, END, freq=bad_freq)
