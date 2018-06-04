# coding=utf-8
# pylint: disable-msg=E1101,W0612

import locale
import calendar
import pytest

from datetime import datetime, date

import numpy as np
import pandas as pd

from pandas.core.dtypes.common import is_integer_dtype, is_list_like
from pandas import (Index, Series, DataFrame, bdate_range,
                    date_range, period_range, timedelta_range,
                    PeriodIndex, DatetimeIndex, TimedeltaIndex)
import pandas.core.common as com

from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesDatetimeValues(TestData):

    def test_dt_namespace_accessor(self):

        # GH 7207, 11128
        # test .dt namespace accessor

        ok_for_period = PeriodIndex._datetimelike_ops
        ok_for_period_methods = ['strftime', 'to_timestamp', 'asfreq']
        ok_for_dt = DatetimeIndex._datetimelike_ops
        ok_for_dt_methods = ['to_period', 'to_pydatetime', 'tz_localize',
                             'tz_convert', 'normalize', 'strftime', 'round',
                             'floor', 'ceil', 'day_name', 'month_name']
        ok_for_td = TimedeltaIndex._datetimelike_ops
        ok_for_td_methods = ['components', 'to_pytimedelta', 'total_seconds',
                             'round', 'floor', 'ceil']

        def get_expected(s, name):
            result = getattr(Index(s._values), prop)
            if isinstance(result, np.ndarray):
                if is_integer_dtype(result):
                    result = result.astype('int64')
            elif not is_list_like(result):
                return result
            return Series(result, index=s.index, name=s.name)

        def compare(s, name):
            a = getattr(s.dt, prop)
            b = get_expected(s, prop)
            if not (is_list_like(a) and is_list_like(b)):
                assert a == b
            else:
                tm.assert_series_equal(a, b)

        # datetimeindex
        cases = [Series(date_range('20130101', periods=5), name='xxx'),
                 Series(date_range('20130101', periods=5, freq='s'),
                        name='xxx'),
                 Series(date_range('20130101 00:00:00', periods=5, freq='ms'),
                        name='xxx')]
        for s in cases:
            for prop in ok_for_dt:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_dt_methods:
                getattr(s.dt, prop)

            result = s.dt.to_pydatetime()
            assert isinstance(result, np.ndarray)
            assert result.dtype == object

            result = s.dt.tz_localize('US/Eastern')
            exp_values = DatetimeIndex(s.values).tz_localize('US/Eastern')
            expected = Series(exp_values, index=s.index, name='xxx')
            tm.assert_series_equal(result, expected)

            tz_result = result.dt.tz
            assert str(tz_result) == 'US/Eastern'
            freq_result = s.dt.freq
            assert freq_result == DatetimeIndex(s.values, freq='infer').freq

            # let's localize, then convert
            result = s.dt.tz_localize('UTC').dt.tz_convert('US/Eastern')
            exp_values = (DatetimeIndex(s.values).tz_localize('UTC')
                                                 .tz_convert('US/Eastern'))
            expected = Series(exp_values, index=s.index, name='xxx')
            tm.assert_series_equal(result, expected)

        # round
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.round('D')
        expected = Series(pd.to_datetime(['2012-01-02', '2012-01-02',
                                          '2012-01-01']), name='xxx')
        tm.assert_series_equal(result, expected)

        # round with tz
        result = (s.dt.tz_localize('UTC')
                   .dt.tz_convert('US/Eastern')
                   .dt.round('D'))
        exp_values = pd.to_datetime(['2012-01-01', '2012-01-01',
                                     '2012-01-01']).tz_localize('US/Eastern')
        expected = Series(exp_values, name='xxx')
        tm.assert_series_equal(result, expected)

        # floor
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.floor('D')
        expected = Series(pd.to_datetime(['2012-01-01', '2012-01-01',
                                          '2012-01-01']), name='xxx')
        tm.assert_series_equal(result, expected)

        # ceil
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.ceil('D')
        expected = Series(pd.to_datetime(['2012-01-02', '2012-01-02',
                                          '2012-01-02']), name='xxx')
        tm.assert_series_equal(result, expected)

        # datetimeindex with tz
        s = Series(date_range('20130101', periods=5, tz='US/Eastern'),
                   name='xxx')
        for prop in ok_for_dt:

            # we test freq below
            if prop != 'freq':
                compare(s, prop)

        for prop in ok_for_dt_methods:
            getattr(s.dt, prop)

        result = s.dt.to_pydatetime()
        assert isinstance(result, np.ndarray)
        assert result.dtype == object

        result = s.dt.tz_convert('CET')
        expected = Series(s._values.tz_convert('CET'),
                          index=s.index, name='xxx')
        tm.assert_series_equal(result, expected)

        tz_result = result.dt.tz
        assert str(tz_result) == 'CET'
        freq_result = s.dt.freq
        assert freq_result == DatetimeIndex(s.values, freq='infer').freq

        # timedelta index
        cases = [Series(timedelta_range('1 day', periods=5),
                        index=list('abcde'), name='xxx'),
                 Series(timedelta_range('1 day 01:23:45', periods=5,
                                        freq='s'), name='xxx'),
                 Series(timedelta_range('2 days 01:23:45.012345', periods=5,
                                        freq='ms'), name='xxx')]
        for s in cases:
            for prop in ok_for_td:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_td_methods:
                getattr(s.dt, prop)

            result = s.dt.components
            assert isinstance(result, DataFrame)
            tm.assert_index_equal(result.index, s.index)

            result = s.dt.to_pytimedelta()
            assert isinstance(result, np.ndarray)
            assert result.dtype == object

            result = s.dt.total_seconds()
            assert isinstance(result, pd.Series)
            assert result.dtype == 'float64'

            freq_result = s.dt.freq
            assert freq_result == TimedeltaIndex(s.values, freq='infer').freq

        # both
        index = date_range('20130101', periods=3, freq='D')
        s = Series(date_range('20140204', periods=3, freq='s'),
                   index=index, name='xxx')
        exp = Series(np.array([2014, 2014, 2014], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.year, exp)

        exp = Series(np.array([2, 2, 2], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.month, exp)

        exp = Series(np.array([0, 1, 2], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.second, exp)

        exp = pd.Series([s[0]] * 3, index=index, name='xxx')
        tm.assert_series_equal(s.dt.normalize(), exp)

        # periodindex
        cases = [Series(period_range('20130101', periods=5, freq='D'),
                        name='xxx')]
        for s in cases:
            for prop in ok_for_period:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_period_methods:
                getattr(s.dt, prop)

            freq_result = s.dt.freq
            assert freq_result == PeriodIndex(s.values).freq

        # test limited display api
        def get_dir(s):
            results = [r for r in s.dt.__dir__() if not r.startswith('_')]
            return list(sorted(set(results)))

        s = Series(date_range('20130101', periods=5, freq='D'), name='xxx')
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_dt + ok_for_dt_methods))))

        s = Series(period_range('20130101', periods=5,
                                freq='D', name='xxx').astype(object))
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_period + ok_for_period_methods))))

        # 11295
        # ambiguous time error on the conversions
        s = Series(pd.date_range('2015-01-01', '2016-01-01',
                                 freq='T'), name='xxx')
        s = s.dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_dt + ok_for_dt_methods))))
        exp_values = pd.date_range('2015-01-01', '2016-01-01', freq='T',
                                   tz='UTC').tz_convert('America/Chicago')
        expected = Series(exp_values, name='xxx')
        tm.assert_series_equal(s, expected)

        # no setting allowed
        s = Series(date_range('20130101', periods=5, freq='D'), name='xxx')
        with tm.assert_raises_regex(ValueError, "modifications"):
            s.dt.hour = 5

        # trying to set a copy
        with pd.option_context('chained_assignment', 'raise'):

            def f():
                s.dt.hour[0] = 5

            pytest.raises(com.SettingWithCopyError, f)

    def test_dt_namespace_accessor_categorical(self):
        # GH 19468
        dti = DatetimeIndex(['20171111', '20181212']).repeat(2)
        s = Series(pd.Categorical(dti), name='foo')
        result = s.dt.year
        expected = Series([2017, 2017, 2018, 2018], name='foo')
        tm.assert_series_equal(result, expected)

    def test_dt_accessor_no_new_attributes(self):
        # https://github.com/pandas-dev/pandas/issues/10673
        s = Series(date_range('20130101', periods=5, freq='D'))
        with tm.assert_raises_regex(AttributeError,
                                    "You cannot add any new attribute"):
            s.dt.xlabel = "a"

    @pytest.mark.parametrize('time_locale', [
        None] if tm.get_locales() is None else [None] + tm.get_locales())
    def test_dt_accessor_datetime_name_accessors(self, time_locale):
        # Test Monday -> Sunday and January -> December, in that sequence
        if time_locale is None:
            # If the time_locale is None, day-name and month_name should
            # return the english attributes
            expected_days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday',
                             'Friday', 'Saturday', 'Sunday']
            expected_months = ['January', 'February', 'March', 'April', 'May',
                               'June', 'July', 'August', 'September',
                               'October', 'November', 'December']
        else:
            with tm.set_locale(time_locale, locale.LC_TIME):
                expected_days = calendar.day_name[:]
                expected_months = calendar.month_name[1:]

        s = Series(DatetimeIndex(freq='D', start=datetime(1998, 1, 1),
                                 periods=365))
        english_days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday',
                        'Friday', 'Saturday', 'Sunday']
        for day, name, eng_name in zip(range(4, 11),
                                       expected_days,
                                       english_days):
            name = name.capitalize()
            assert s.dt.weekday_name[day] == eng_name
            assert s.dt.day_name(locale=time_locale)[day] == name
        s = s.append(Series([pd.NaT]))
        assert np.isnan(s.dt.day_name(locale=time_locale).iloc[-1])

        s = Series(DatetimeIndex(freq='M', start='2012', end='2013'))
        result = s.dt.month_name(locale=time_locale)
        expected = Series([month.capitalize() for month in expected_months])
        tm.assert_series_equal(result, expected)
        for s_date, expected in zip(s, expected_months):
            result = s_date.month_name(locale=time_locale)
            assert result == expected.capitalize()
        s = s.append(Series([pd.NaT]))
        assert np.isnan(s.dt.month_name(locale=time_locale).iloc[-1])

    def test_strftime(self):
        # GH 10086
        s = Series(date_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03',
                           '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('2015-02-03 11:22:33.4567', periods=5))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/04 11-22-33',
                           '2015/02/05 11-22-33', '2015/02/06 11-22-33',
                           '2015/02/07 11-22-33'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03',
                           '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range(
            '2015-02-03 11:22:33.4567', periods=5, freq='s'))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/03 11-22-34',
                           '2015/02/03 11-22-35', '2015/02/03 11-22-36',
                           '2015/02/03 11-22-37'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('20130101', periods=5))
        s.iloc[0] = pd.NaT
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['NaT', '2013/01/02', '2013/01/03', '2013/01/04',
                           '2013/01/05'])
        tm.assert_series_equal(result, expected)

        datetime_index = date_range('20150301', periods=5)
        result = datetime_index.strftime("%Y/%m/%d")

        expected = Index(['2015/03/01', '2015/03/02', '2015/03/03',
                          '2015/03/04', '2015/03/05'], dtype=np.object_)
        # dtype may be S10 or U10 depending on python version
        tm.assert_index_equal(result, expected)

        period_index = period_range('20150301', periods=5)
        result = period_index.strftime("%Y/%m/%d")
        expected = Index(['2015/03/01', '2015/03/02', '2015/03/03',
                          '2015/03/04', '2015/03/05'], dtype='=U10')
        tm.assert_index_equal(result, expected)

        s = Series([datetime(2013, 1, 1, 2, 32, 59), datetime(2013, 1, 2, 14,
                                                              32, 1)])
        result = s.dt.strftime('%Y-%m-%d %H:%M:%S')
        expected = Series(["2013-01-01 02:32:59", "2013-01-02 14:32:01"])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=4, freq='H'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S')
        expected = Series(["2013/01/01 00:00:00", "2013/01/01 01:00:00",
                           "2013/01/01 02:00:00", "2013/01/01 03:00:00"])

        s = Series(period_range('20130101', periods=4, freq='L'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S.%l')
        expected = Series(["2013/01/01 00:00:00.000",
                           "2013/01/01 00:00:00.001",
                           "2013/01/01 00:00:00.002",
                           "2013/01/01 00:00:00.003"])
        tm.assert_series_equal(result, expected)

    def test_valid_dt_with_missing_values(self):

        from datetime import date, time

        # GH 8689
        s = Series(date_range('20130101', periods=5, freq='D'))
        s.iloc[2] = pd.NaT

        for attr in ['microsecond', 'nanosecond', 'second', 'minute', 'hour',
                     'day']:
            expected = getattr(s.dt, attr).copy()
            expected.iloc[2] = np.nan
            result = getattr(s.dt, attr)
            tm.assert_series_equal(result, expected)

        result = s.dt.date
        expected = Series(
            [date(2013, 1, 1), date(2013, 1, 2), np.nan, date(2013, 1, 4),
             date(2013, 1, 5)], dtype='object')
        tm.assert_series_equal(result, expected)

        result = s.dt.time
        expected = Series(
            [time(0), time(0), np.nan, time(0), time(0)], dtype='object')
        tm.assert_series_equal(result, expected)

    def test_dt_accessor_api(self):
        # GH 9322
        from pandas.core.indexes.accessors import (
            CombinedDatetimelikeProperties, DatetimeProperties)
        assert Series.dt is CombinedDatetimelikeProperties

        s = Series(date_range('2000-01-01', periods=3))
        assert isinstance(s.dt, DatetimeProperties)

        for s in [Series(np.arange(5)), Series(list('abcde')),
                  Series(np.random.randn(5))]:
            with tm.assert_raises_regex(AttributeError,
                                        "only use .dt accessor"):
                s.dt
            assert not hasattr(s, 'dt')

    def test_between(self):
        s = Series(bdate_range('1/1/2000', periods=20).astype(object))
        s[::2] = np.nan

        result = s[s.between(s[3], s[17])]
        expected = s[3:18].dropna()
        assert_series_equal(result, expected)

        result = s[s.between(s[3], s[17], inclusive=False)]
        expected = s[5:16].dropna()
        assert_series_equal(result, expected)

    def test_date_tz(self):
        # GH11757
        rng = pd.DatetimeIndex(['2014-04-04 23:56',
                                '2014-07-18 21:24',
                                '2015-11-22 22:14'], tz="US/Eastern")
        s = Series(rng)
        expected = Series([date(2014, 4, 4),
                           date(2014, 7, 18),
                           date(2015, 11, 22)])
        assert_series_equal(s.dt.date, expected)
        assert_series_equal(s.apply(lambda x: x.date()), expected)

    def test_datetime_understood(self):
        # Ensures it doesn't fail to create the right series
        # reported in issue#16726
        series = pd.Series(pd.date_range("2012-01-01", periods=3))
        offset = pd.offsets.DateOffset(days=6)
        result = series - offset
        expected = pd.Series(pd.to_datetime([
            '2011-12-26', '2011-12-27', '2011-12-28']))
        tm.assert_series_equal(result, expected)
