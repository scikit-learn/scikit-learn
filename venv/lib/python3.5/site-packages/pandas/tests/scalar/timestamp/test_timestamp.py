""" test the scalar Timestamp """

import pytz
import pytest
import dateutil
import calendar
import locale
import numpy as np

from dateutil.tz import tzutc
from pytz import timezone, utc
from datetime import datetime, timedelta

import pandas.util.testing as tm
import pandas.util._test_decorators as td

from pandas.tseries import offsets

from pandas._libs.tslibs import conversion
from pandas._libs.tslibs.timezones import get_timezone, dateutil_gettz as gettz

from pandas.errors import OutOfBoundsDatetime
from pandas.compat import long, PY3
from pandas.compat.numpy import np_datetime64_compat
from pandas import Timestamp, Period, Timedelta, NaT


class TestTimestampProperties(object):

    def test_properties_business(self):
        ts = Timestamp('2017-10-01', freq='B')
        control = Timestamp('2017-10-01')
        assert ts.dayofweek == 6
        assert not ts.is_month_start    # not a weekday
        assert not ts.is_quarter_start  # not a weekday
        # Control case: non-business is month/qtr start
        assert control.is_month_start
        assert control.is_quarter_start

        ts = Timestamp('2017-09-30', freq='B')
        control = Timestamp('2017-09-30')
        assert ts.dayofweek == 5
        assert not ts.is_month_end    # not a weekday
        assert not ts.is_quarter_end  # not a weekday
        # Control case: non-business is month/qtr start
        assert control.is_month_end
        assert control.is_quarter_end

    def test_fields(self):
        def check(value, equal):
            # that we are int/long like
            assert isinstance(value, (int, long))
            assert value == equal

        # GH 10050
        ts = Timestamp('2015-05-10 09:06:03.000100001')
        check(ts.year, 2015)
        check(ts.month, 5)
        check(ts.day, 10)
        check(ts.hour, 9)
        check(ts.minute, 6)
        check(ts.second, 3)
        pytest.raises(AttributeError, lambda: ts.millisecond)
        check(ts.microsecond, 100)
        check(ts.nanosecond, 1)
        check(ts.dayofweek, 6)
        check(ts.quarter, 2)
        check(ts.dayofyear, 130)
        check(ts.week, 19)
        check(ts.daysinmonth, 31)
        check(ts.daysinmonth, 31)

        # GH 13303
        ts = Timestamp('2014-12-31 23:59:00-05:00', tz='US/Eastern')
        check(ts.year, 2014)
        check(ts.month, 12)
        check(ts.day, 31)
        check(ts.hour, 23)
        check(ts.minute, 59)
        check(ts.second, 0)
        pytest.raises(AttributeError, lambda: ts.millisecond)
        check(ts.microsecond, 0)
        check(ts.nanosecond, 0)
        check(ts.dayofweek, 2)
        check(ts.quarter, 4)
        check(ts.dayofyear, 365)
        check(ts.week, 1)
        check(ts.daysinmonth, 31)

        ts = Timestamp('2014-01-01 00:00:00+01:00')
        starts = ['is_month_start', 'is_quarter_start', 'is_year_start']
        for start in starts:
            assert getattr(ts, start)
        ts = Timestamp('2014-12-31 23:59:59+01:00')
        ends = ['is_month_end', 'is_year_end', 'is_quarter_end']
        for end in ends:
            assert getattr(ts, end)

    # GH 12806
    @pytest.mark.parametrize('data',
                             [Timestamp('2017-08-28 23:00:00'),
                              Timestamp('2017-08-28 23:00:00', tz='EST')])
    @pytest.mark.parametrize('time_locale', [
        None] if tm.get_locales() is None else [None] + tm.get_locales())
    def test_names(self, data, time_locale):
        # GH 17354
        # Test .weekday_name, .day_name(), .month_name
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            assert data.weekday_name == 'Monday'
        if time_locale is None:
            expected_day = 'Monday'
            expected_month = 'August'
        else:
            with tm.set_locale(time_locale, locale.LC_TIME):
                expected_day = calendar.day_name[0].capitalize()
                expected_month = calendar.month_name[8].capitalize()

        assert data.day_name(time_locale) == expected_day
        assert data.month_name(time_locale) == expected_month

        # Test NaT
        nan_ts = Timestamp(NaT)
        assert np.isnan(nan_ts.day_name(time_locale))
        assert np.isnan(nan_ts.month_name(time_locale))

    def test_is_leap_year(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH 13727
        dt = Timestamp('2000-01-01 00:00:00', tz=tz)
        assert dt.is_leap_year
        assert isinstance(dt.is_leap_year, bool)

        dt = Timestamp('1999-01-01 00:00:00', tz=tz)
        assert not dt.is_leap_year

        dt = Timestamp('2004-01-01 00:00:00', tz=tz)
        assert dt.is_leap_year

        dt = Timestamp('2100-01-01 00:00:00', tz=tz)
        assert not dt.is_leap_year

    def test_woy_boundary(self):
        # make sure weeks at year boundaries are correct
        d = datetime(2013, 12, 31)
        result = Timestamp(d).week
        expected = 1  # ISO standard
        assert result == expected

        d = datetime(2008, 12, 28)
        result = Timestamp(d).week
        expected = 52  # ISO standard
        assert result == expected

        d = datetime(2009, 12, 31)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        d = datetime(2010, 1, 1)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        d = datetime(2010, 1, 3)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        result = np.array([Timestamp(datetime(*args)).week
                           for args in [(2000, 1, 1), (2000, 1, 2), (
                               2005, 1, 1), (2005, 1, 2)]])
        assert (result == [52, 52, 53, 53]).all()


class TestTimestampConstructors(object):

    def test_constructor(self):
        base_str = '2014-07-01 09:00'
        base_dt = datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        assert (calendar.timegm(base_dt.timetuple()) * 1000000000 ==
                base_expected)

        tests = [(base_str, base_dt, base_expected),
                 ('2014-07-01 10:00', datetime(2014, 7, 1, 10),
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 09:00:00.000008000',
                  datetime(2014, 7, 1, 9, 0, 0, 8),
                  base_expected + 8000),
                 ('2014-07-01 09:00:00.000000005',
                  Timestamp('2014-07-01 09:00:00.000000005'),
                  base_expected + 5)]

        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0), ('Asia/Tokyo', 9),
                     ('US/Eastern', -4), ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3),
                     (dateutil.tz.tzoffset(None, 18000), 5)]

        for date_str, date, expected in tests:
            for result in [Timestamp(date_str), Timestamp(date)]:
                # only with timestring
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

            # with timezone
            for tz, offset in timezones:
                for result in [Timestamp(date_str, tz=tz), Timestamp(date,
                                                                     tz=tz)]:
                    expected_tz = expected - offset * 3600 * 1000000000
                    assert result.value == expected_tz
                    assert conversion.pydt_to_i8(result) == expected_tz

                    # should preserve tz
                    result = Timestamp(result)
                    assert result.value == expected_tz
                    assert conversion.pydt_to_i8(result) == expected_tz

                    # should convert to UTC
                    result = Timestamp(result, tz='UTC')
                    expected_utc = expected - offset * 3600 * 1000000000
                    assert result.value == expected_utc
                    assert conversion.pydt_to_i8(result) == expected_utc

    def test_constructor_with_stringoffset(self):
        # GH 7833
        base_str = '2014-07-01 11:00:00+02:00'
        base_dt = datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        assert (calendar.timegm(base_dt.timetuple()) * 1000000000 ==
                base_expected)

        tests = [(base_str, base_expected),
                 ('2014-07-01 12:00:00+02:00',
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 11:00:00.000008000+02:00', base_expected + 8000),
                 ('2014-07-01 11:00:00.000000005+02:00', base_expected + 5)]

        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0), ('Asia/Tokyo', 9),
                     ('US/Eastern', -4), ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3),
                     (dateutil.tz.tzoffset(None, 18000), 5)]

        for date_str, expected in tests:
            for result in [Timestamp(date_str)]:
                # only with timestring
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

            # with timezone
            for tz, offset in timezones:
                result = Timestamp(date_str, tz=tz)
                expected_tz = expected
                assert result.value == expected_tz
                assert conversion.pydt_to_i8(result) == expected_tz

                # should preserve tz
                result = Timestamp(result)
                assert result.value == expected_tz
                assert conversion.pydt_to_i8(result) == expected_tz

                # should convert to UTC
                result = Timestamp(result, tz='UTC')
                expected_utc = expected
                assert result.value == expected_utc
                assert conversion.pydt_to_i8(result) == expected_utc

        # This should be 2013-11-01 05:00 in UTC
        # converted to Chicago tz
        result = Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')
        assert result.value == Timestamp('2013-11-01 05:00').value
        expected = "Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')"  # noqa
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2013-11-01 05:00 in UTC
        # converted to Tokyo tz (+09:00)
        result = Timestamp('2013-11-01 00:00:00-0500', tz='Asia/Tokyo')
        assert result.value == Timestamp('2013-11-01 05:00').value
        expected = "Timestamp('2013-11-01 14:00:00+0900', tz='Asia/Tokyo')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # GH11708
        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Katmandu
        result = Timestamp("2015-11-18 15:45:00+05:45", tz="Asia/Katmandu")
        assert result.value == Timestamp("2015-11-18 10:00").value
        expected = "Timestamp('2015-11-18 15:45:00+0545', tz='Asia/Katmandu')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Kolkata
        result = Timestamp("2015-11-18 15:30:00+05:30", tz="Asia/Kolkata")
        assert result.value == Timestamp("2015-11-18 10:00").value
        expected = "Timestamp('2015-11-18 15:30:00+0530', tz='Asia/Kolkata')"
        assert repr(result) == expected
        assert result == eval(repr(result))

    def test_constructor_invalid(self):
        with tm.assert_raises_regex(TypeError, 'Cannot convert input'):
            Timestamp(slice(2))
        with tm.assert_raises_regex(ValueError, 'Cannot convert Period'):
            Timestamp(Period('1000-01-01'))

    def test_constructor_invalid_tz(self):
        # GH#17690
        with tm.assert_raises_regex(TypeError, 'must be a datetime.tzinfo'):
            Timestamp('2017-10-22', tzinfo='US/Eastern')

        with tm.assert_raises_regex(ValueError, 'at most one of'):
            Timestamp('2017-10-22', tzinfo=utc, tz='UTC')

        with tm.assert_raises_regex(ValueError, "Invalid frequency:"):
            # GH#5168
            # case where user tries to pass tz as an arg, not kwarg, gets
            # interpreted as a `freq`
            Timestamp('2012-01-01', 'US/Pacific')

    def test_constructor_tz_or_tzinfo(self):
        # GH#17943, GH#17690, GH#5168
        stamps = [Timestamp(year=2017, month=10, day=22, tz='UTC'),
                  Timestamp(year=2017, month=10, day=22, tzinfo=utc),
                  Timestamp(year=2017, month=10, day=22, tz=utc),
                  Timestamp(datetime(2017, 10, 22), tzinfo=utc),
                  Timestamp(datetime(2017, 10, 22), tz='UTC'),
                  Timestamp(datetime(2017, 10, 22), tz=utc)]
        assert all(ts == stamps[0] for ts in stamps)

    def test_constructor_positional(self):
        # see gh-10758
        with pytest.raises(TypeError):
            Timestamp(2000, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 0, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 13, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 1, 0)
        with pytest.raises(ValueError):
            Timestamp(2000, 1, 32)

        # see gh-11630
        assert (repr(Timestamp(2015, 11, 12)) ==
                repr(Timestamp('20151112')))
        assert (repr(Timestamp(2015, 11, 12, 1, 2, 3, 999999)) ==
                repr(Timestamp('2015-11-12 01:02:03.999999')))

    def test_constructor_keyword(self):
        # GH 10758
        with pytest.raises(TypeError):
            Timestamp(year=2000, month=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=0, day=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=13, day=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=1, day=0)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=1, day=32)

        assert (repr(Timestamp(year=2015, month=11, day=12)) ==
                repr(Timestamp('20151112')))

        assert (repr(Timestamp(year=2015, month=11, day=12, hour=1, minute=2,
                               second=3, microsecond=999999)) ==
                repr(Timestamp('2015-11-12 01:02:03.999999')))

    def test_constructor_fromordinal(self):
        base = datetime(2000, 1, 1)

        ts = Timestamp.fromordinal(base.toordinal(), freq='D')
        assert base == ts
        assert ts.freq == 'D'
        assert base.toordinal() == ts.toordinal()

        ts = Timestamp.fromordinal(base.toordinal(), tz='US/Eastern')
        assert Timestamp('2000-01-01', tz='US/Eastern') == ts
        assert base.toordinal() == ts.toordinal()

        # GH#3042
        dt = datetime(2011, 4, 16, 0, 0)
        ts = Timestamp.fromordinal(dt.toordinal())
        assert ts.to_pydatetime() == dt

        # with a tzinfo
        stamp = Timestamp('2011-4-16', tz='US/Eastern')
        dt_tz = stamp.to_pydatetime()
        ts = Timestamp.fromordinal(dt_tz.toordinal(), tz='US/Eastern')
        assert ts.to_pydatetime() == dt_tz

    @pytest.mark.parametrize('result', [
        Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6), nanosecond=1),
        Timestamp(year=2000, month=1, day=2, hour=3, minute=4, second=5,
                  microsecond=6, nanosecond=1),
        Timestamp(year=2000, month=1, day=2, hour=3, minute=4, second=5,
                  microsecond=6, nanosecond=1, tz='UTC'),
        Timestamp(2000, 1, 2, 3, 4, 5, 6, 1, None),
        Timestamp(2000, 1, 2, 3, 4, 5, 6, 1, pytz.UTC)])
    def test_constructor_nanosecond(self, result):
        # GH 18898
        expected = Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6), tz=result.tz)
        expected = expected + Timedelta(nanoseconds=1)
        assert result == expected

    @pytest.mark.parametrize('arg', ['year', 'month', 'day', 'hour', 'minute',
                                     'second', 'microsecond', 'nanosecond'])
    def test_invalid_date_kwarg_with_string_input(self, arg):
        kwarg = {arg: 1}
        with pytest.raises(ValueError):
            Timestamp('2010-10-10 12:59:59.999999999', **kwarg)

    def test_out_of_bounds_value(self):
        one_us = np.timedelta64(1).astype('timedelta64[us]')

        # By definition we can't go out of bounds in [ns], so we
        # convert the datetime64s to [us] so we can go out of bounds
        min_ts_us = np.datetime64(Timestamp.min).astype('M8[us]')
        max_ts_us = np.datetime64(Timestamp.max).astype('M8[us]')

        # No error for the min/max datetimes
        Timestamp(min_ts_us)
        Timestamp(max_ts_us)

        # One us less than the minimum is an error
        with pytest.raises(ValueError):
            Timestamp(min_ts_us - one_us)

        # One us more than the maximum is an error
        with pytest.raises(ValueError):
            Timestamp(max_ts_us + one_us)

    def test_out_of_bounds_string(self):
        with pytest.raises(ValueError):
            Timestamp('1676-01-01')
        with pytest.raises(ValueError):
            Timestamp('2263-01-01')

    def test_barely_out_of_bounds(self):
        # GH#19529
        # GH#19382 close enough to bounds that dropping nanos would result
        # in an in-bounds datetime
        with pytest.raises(OutOfBoundsDatetime):
            Timestamp('2262-04-11 23:47:16.854775808')

    def test_bounds_with_different_units(self):
        out_of_bounds_dates = ('1677-09-21', '2262-04-12')

        time_units = ('D', 'h', 'm', 's', 'ms', 'us')

        for date_string in out_of_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, dtype='M8[%s]' % unit)
                with pytest.raises(ValueError):
                    Timestamp(dt64)

        in_bounds_dates = ('1677-09-23', '2262-04-11')

        for date_string in in_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, dtype='M8[%s]' % unit)
                Timestamp(dt64)

    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        Timestamp(Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        Timestamp(Timestamp.max)

    def test_now(self):
        # GH#9000
        ts_from_string = Timestamp('now')
        ts_from_method = Timestamp.now()
        ts_datetime = datetime.now()

        ts_from_string_tz = Timestamp('now', tz='US/Eastern')
        ts_from_method_tz = Timestamp.now(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        assert abs(ts_from_method - ts_from_string) < delta
        assert abs(ts_datetime - ts_from_method) < delta
        assert abs(ts_from_method_tz - ts_from_string_tz) < delta
        assert (abs(ts_from_string_tz.tz_localize(None) -
                    ts_from_method_tz.tz_localize(None)) < delta)

    def test_today(self):
        ts_from_string = Timestamp('today')
        ts_from_method = Timestamp.today()
        ts_datetime = datetime.today()

        ts_from_string_tz = Timestamp('today', tz='US/Eastern')
        ts_from_method_tz = Timestamp.today(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        assert abs(ts_from_method - ts_from_string) < delta
        assert abs(ts_datetime - ts_from_method) < delta
        assert abs(ts_from_method_tz - ts_from_string_tz) < delta
        assert (abs(ts_from_string_tz.tz_localize(None) -
                    ts_from_method_tz.tz_localize(None)) < delta)

    @pytest.mark.parametrize('tz', [None, pytz.timezone('US/Pacific')])
    def test_disallow_setting_tz(self, tz):
        # GH 3746
        ts = Timestamp('2010')
        with pytest.raises(AttributeError):
            ts.tz = tz


class TestTimestamp(object):

    def test_tz(self):
        tstr = '2014-02-01 09:00'
        ts = Timestamp(tstr)
        local = ts.tz_localize('Asia/Tokyo')
        assert local.hour == 9
        assert local == Timestamp(tstr, tz='Asia/Tokyo')
        conv = local.tz_convert('US/Eastern')
        assert conv == Timestamp('2014-01-31 19:00', tz='US/Eastern')
        assert conv.hour == 19

        # preserves nanosecond
        ts = Timestamp(tstr) + offsets.Nano(5)
        local = ts.tz_localize('Asia/Tokyo')
        assert local.hour == 9
        assert local.nanosecond == 5
        conv = local.tz_convert('US/Eastern')
        assert conv.nanosecond == 5
        assert conv.hour == 19

    def test_utc_z_designator(self):
        assert get_timezone(Timestamp('2014-11-02 01:00Z').tzinfo) == 'UTC'

    def test_asm8(self):
        np.random.seed(7960929)
        ns = [Timestamp.min.value, Timestamp.max.value, 1000]

        for n in ns:
            assert (Timestamp(n).asm8.view('i8') ==
                    np.datetime64(n, 'ns').view('i8') == n)

        assert (Timestamp('nat').asm8.view('i8') ==
                np.datetime64('nat', 'ns').view('i8'))

    def test_class_ops_pytz(self):
        def compare(x, y):
            assert (int(Timestamp(x).value / 1e9) ==
                    int(Timestamp(y).value / 1e9))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(timezone('UTC')))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_class_ops_dateutil(self):
        def compare(x, y):
            assert (int(np.round(Timestamp(x).value / 1e9)) ==
                    int(np.round(Timestamp(y).value / 1e9)))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(tzutc()))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_basics_nanos(self):
        val = np.int64(946684800000000000).view('M8[ns]')
        stamp = Timestamp(val.view('i8') + 500)
        assert stamp.year == 2000
        assert stamp.month == 1
        assert stamp.microsecond == 0
        assert stamp.nanosecond == 500

        # GH 14415
        val = np.iinfo(np.int64).min + 80000000000000
        stamp = Timestamp(val)
        assert stamp.year == 1677
        assert stamp.month == 9
        assert stamp.day == 21
        assert stamp.microsecond == 145224
        assert stamp.nanosecond == 192

    def test_unit(self):

        def check(val, unit=None, h=1, s=1, us=0):
            stamp = Timestamp(val, unit=unit)
            assert stamp.year == 2000
            assert stamp.month == 1
            assert stamp.day == 1
            assert stamp.hour == h
            if unit != 'D':
                assert stamp.minute == 1
                assert stamp.second == s
                assert stamp.microsecond == us
            else:
                assert stamp.minute == 0
                assert stamp.second == 0
                assert stamp.microsecond == 0
            assert stamp.nanosecond == 0

        ts = Timestamp('20000101 01:01:01')
        val = ts.value
        days = (ts - Timestamp('1970-01-01')).days

        check(val)
        check(val / long(1000), unit='us')
        check(val / long(1000000), unit='ms')
        check(val / long(1000000000), unit='s')
        check(days, unit='D', h=0)

        # using truediv, so these are like floats
        if PY3:
            check((val + 500000) / long(1000000000), unit='s', us=500)
            check((val + 500000000) / long(1000000000), unit='s', us=500000)
            check((val + 500000) / long(1000000), unit='ms', us=500)

        # get chopped in py2
        else:
            check((val + 500000) / long(1000000000), unit='s')
            check((val + 500000000) / long(1000000000), unit='s')
            check((val + 500000) / long(1000000), unit='ms')

        # ok
        check((val + 500000) / long(1000), unit='us', us=500)
        check((val + 500000000) / long(1000000), unit='ms', us=500000)

        # floats
        check(val / 1000.0 + 5, unit='us', us=5)
        check(val / 1000.0 + 5000, unit='us', us=5000)
        check(val / 1000000.0 + 0.5, unit='ms', us=500)
        check(val / 1000000.0 + 0.005, unit='ms', us=5)
        check(val / 1000000000.0 + 0.5, unit='s', us=500000)
        check(days + 0.5, unit='D', h=12)

    def test_roundtrip(self):

        # test value to string and back conversions
        # further test accessors
        base = Timestamp('20140101 00:00:00')

        result = Timestamp(base.value + Timedelta('5ms').value)
        assert result == Timestamp(str(base) + ".005000")
        assert result.microsecond == 5000

        result = Timestamp(base.value + Timedelta('5us').value)
        assert result == Timestamp(str(base) + ".000005")
        assert result.microsecond == 5

        result = Timestamp(base.value + Timedelta('5ns').value)
        assert result == Timestamp(str(base) + ".000000005")
        assert result.nanosecond == 5
        assert result.microsecond == 0

        result = Timestamp(base.value + Timedelta('6ms 5us').value)
        assert result == Timestamp(str(base) + ".006005")
        assert result.microsecond == 5 + 6 * 1000

        result = Timestamp(base.value + Timedelta('200ms 5us').value)
        assert result == Timestamp(str(base) + ".200005")
        assert result.microsecond == 5 + 200 * 1000

    def test_hash_equivalent(self):
        d = {datetime(2011, 1, 1): 5}
        stamp = Timestamp(datetime(2011, 1, 1))
        assert d[stamp] == 5


class TestTimestampNsOperations(object):

    def setup_method(self, method):
        self.timestamp = Timestamp(datetime.utcnow())

    def assert_ns_timedelta(self, modified_timestamp, expected_value):
        value = self.timestamp.value
        modified_value = modified_timestamp.value

        assert modified_value - value == expected_value

    def test_timedelta_ns_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ns'),
                                 -123)

    def test_timedelta_ns_based_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(
            1234567898, 'ns'), 1234567898)

    def test_timedelta_us_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'us'),
                                 -123000)

    def test_timedelta_ms_arithmetic(self):
        time = self.timestamp + np.timedelta64(-123, 'ms')
        self.assert_ns_timedelta(time, -123000000)

    def test_nanosecond_string_parsing(self):
        ts = Timestamp('2013-05-01 07:15:45.123456789')
        # GH 7878
        expected_repr = '2013-05-01 07:15:45.123456789'
        expected_value = 1367392545123456789
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

        ts = Timestamp('2013-05-01 07:15:45.123456789+09:00', tz='Asia/Tokyo')
        assert ts.value == expected_value - 9 * 3600 * 1000000000
        assert expected_repr in repr(ts)

        ts = Timestamp('2013-05-01 07:15:45.123456789', tz='UTC')
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

        ts = Timestamp('2013-05-01 07:15:45.123456789', tz='US/Eastern')
        assert ts.value == expected_value + 4 * 3600 * 1000000000
        assert expected_repr in repr(ts)

        # GH 10041
        ts = Timestamp('20130501T071545.123456789')
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

    def test_nanosecond_timestamp(self):
        # GH 7610
        expected = 1293840000000000005
        t = Timestamp('2011-01-01') + offsets.Nano(5)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        t = Timestamp(t)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        t = Timestamp(np_datetime64_compat('2011-01-01 00:00:00.000000005Z'))
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        expected = 1293840000000000010
        t = t + offsets.Nano(5)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10

        t = Timestamp(t)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10

        t = Timestamp(np_datetime64_compat('2011-01-01 00:00:00.000000010Z'))
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10


class TestTimestampToJulianDate(object):

    def test_compare_1700(self):
        r = Timestamp('1700-06-23').to_julian_date()
        assert r == 2342145.5

    def test_compare_2000(self):
        r = Timestamp('2000-04-12').to_julian_date()
        assert r == 2451646.5

    def test_compare_2100(self):
        r = Timestamp('2100-08-12').to_julian_date()
        assert r == 2488292.5

    def test_compare_hour01(self):
        r = Timestamp('2000-08-12T01:00:00').to_julian_date()
        assert r == 2451768.5416666666666666

    def test_compare_hour13(self):
        r = Timestamp('2000-08-12T13:00:00').to_julian_date()
        assert r == 2451769.0416666666666666


class TestTimestampConversion(object):
    def test_conversion(self):
        # GH#9255
        ts = Timestamp('2000-01-01')

        result = ts.to_pydatetime()
        expected = datetime(2000, 1, 1)
        assert result == expected
        assert type(result) == type(expected)

        result = ts.to_datetime64()
        expected = np.datetime64(ts.value, 'ns')
        assert result == expected
        assert type(result) == type(expected)
        assert result.dtype == expected.dtype

    def test_to_pydatetime_nonzero_nano(self):
        ts = Timestamp('2011-01-01 9:00:00.123456789')

        # Warn the user of data loss (nanoseconds).
        with tm.assert_produces_warning(UserWarning,
                                        check_stacklevel=False):
            expected = datetime(2011, 1, 1, 9, 0, 0, 123456)
            result = ts.to_pydatetime()
            assert result == expected

    def test_timestamp_to_datetime(self):
        stamp = Timestamp('20090415', tz='US/Eastern', freq='D')
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_timestamp_to_datetime_dateutil(self):
        stamp = Timestamp('20090415', tz='dateutil/US/Eastern', freq='D')
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_timestamp_to_datetime_explicit_pytz(self):
        stamp = Timestamp('20090415', tz=pytz.timezone('US/Eastern'), freq='D')
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    @td.skip_if_windows_python_3
    def test_timestamp_to_datetime_explicit_dateutil(self):
        stamp = Timestamp('20090415', tz=gettz('US/Eastern'), freq='D')
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_to_datetime_bijective(self):
        # Ensure that converting to datetime and back only loses precision
        # by going from nanoseconds to microseconds.
        exp_warning = None if Timestamp.max.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            assert (Timestamp(Timestamp.max.to_pydatetime()).value / 1000 ==
                    Timestamp.max.value / 1000)

        exp_warning = None if Timestamp.min.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            assert (Timestamp(Timestamp.min.to_pydatetime()).value / 1000 ==
                    Timestamp.min.value / 1000)
