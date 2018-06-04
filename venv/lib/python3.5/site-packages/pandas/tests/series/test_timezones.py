# -*- coding: utf-8 -*-
"""
Tests for Series timezone-related methods
"""
from datetime import datetime

import pytest
import pytz
import numpy as np
from dateutil.tz import tzoffset

import pandas.util.testing as tm
from pandas._libs import tslib
from pandas._libs.tslibs import timezones
from pandas.compat import lrange
from pandas.core.indexes.datetimes import date_range
from pandas import Series, Timestamp, DatetimeIndex, Index


class TestSeriesTimezones(object):
    # -----------------------------------------------------------------
    # Series.tz_localize
    def test_series_tz_localize(self):

        rng = date_range('1/1/2011', periods=100, freq='H')
        ts = Series(1, index=rng)

        result = ts.tz_localize('utc')
        assert result.index.tz.zone == 'UTC'

        # Can't localize if already tz-aware
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        ts = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, 'Already tz-aware',
                               ts.tz_localize, 'US/Eastern')

    def test_series_tz_localize_ambiguous_bool(self):
        # make sure that we are correctly accepting bool values as ambiguous

        # GH#14402
        ts = Timestamp('2015-11-01 01:00:03')
        expected0 = Timestamp('2015-11-01 01:00:03-0500', tz='US/Central')
        expected1 = Timestamp('2015-11-01 01:00:03-0600', tz='US/Central')

        ser = Series([ts])
        expected0 = Series([expected0])
        expected1 = Series([expected1])

        with pytest.raises(pytz.AmbiguousTimeError):
            ser.dt.tz_localize('US/Central')

        result = ser.dt.tz_localize('US/Central', ambiguous=True)
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize('US/Central', ambiguous=[True])
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize('US/Central', ambiguous=False)
        tm.assert_series_equal(result, expected1)

        result = ser.dt.tz_localize('US/Central', ambiguous=[False])
        tm.assert_series_equal(result, expected1)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_series_tz_localize_empty(self, tzstr):
        # GH#2248
        ser = Series()

        ser2 = ser.tz_localize('utc')
        assert ser2.index.tz == pytz.utc

        ser2 = ser.tz_localize(tzstr)
        timezones.tz_compare(ser2.index.tz, timezones.maybe_get_tz(tzstr))

    # -----------------------------------------------------------------
    # Series.tz_convert

    def test_series_tz_convert(self):
        rng = date_range('1/1/2011', periods=200, freq='D', tz='US/Eastern')
        ts = Series(1, index=rng)

        result = ts.tz_convert('Europe/Berlin')
        assert result.index.tz.zone == 'Europe/Berlin'

        # can't convert tz-naive
        rng = date_range('1/1/2011', periods=200, freq='D')
        ts = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, "Cannot convert tz-naive",
                               ts.tz_convert, 'US/Eastern')

    def test_series_tz_convert_to_utc(self):
        base = DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03'],
                             tz='UTC')
        idx1 = base.tz_convert('Asia/Tokyo')[:2]
        idx2 = base.tz_convert('US/Eastern')[1:]

        res = Series([1, 2], index=idx1) + Series([1, 1], index=idx2)
        tm.assert_series_equal(res, Series([np.nan, 3, np.nan], index=base))

    # -----------------------------------------------------------------
    # Series.append

    def test_series_append_aware(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ser1 = Series([1], index=rng1)
        ser2 = Series([2], index=rng2)
        ts_result = ser1.append(ser2)

        exp_index = DatetimeIndex(['2011-01-01 01:00', '2011-01-01 02:00'],
                                  tz='US/Eastern')
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H', tz='UTC')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H', tz='UTC')
        ser1 = Series([1], index=rng1)
        ser2 = Series([2], index=rng2)
        ts_result = ser1.append(ser2)

        exp_index = DatetimeIndex(['2011-01-01 01:00', '2011-01-01 02:00'],
                                  tz='UTC')
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        utc = rng1.tz
        assert utc == ts_result.index.tz

        # GH#7795
        # different tz coerces to object dtype, not UTC
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Central')
        ser1 = Series([1], index=rng1)
        ser2 = Series([2], index=rng2)
        ts_result = ser1.append(ser2)
        exp_index = Index([Timestamp('1/1/2011 01:00', tz='US/Eastern'),
                           Timestamp('1/1/2011 02:00', tz='US/Central')])
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)

    def test_series_append_aware_naive(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ser1 = Series(np.random.randn(len(rng1)), index=rng1)
        ser2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ser1.append(ser2)

        expected = ser1.index.astype(object).append(ser2.index.astype(object))
        assert ts_result.index.equals(expected)

        # mixed
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = lrange(100)
        ser1 = Series(np.random.randn(len(rng1)), index=rng1)
        ser2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ser1.append(ser2)

        expected = ser1.index.astype(object).append(ser2.index)
        assert ts_result.index.equals(expected)

    def test_series_append_dst(self):
        rng1 = date_range('1/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('8/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        ser1 = Series([1, 2, 3], index=rng1)
        ser2 = Series([10, 11, 12], index=rng2)
        ts_result = ser1.append(ser2)

        exp_index = DatetimeIndex(['2016-01-01 01:00', '2016-01-01 02:00',
                                   '2016-01-01 03:00', '2016-08-01 01:00',
                                   '2016-08-01 02:00', '2016-08-01 03:00'],
                                  tz='US/Eastern')
        exp = Series([1, 2, 3, 10, 11, 12], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

    # -----------------------------------------------------------------

    def test_dateutil_tzoffset_support(self):
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [datetime(2012, 5, 11, 11, tzinfo=tzinfo),
                 datetime(2012, 5, 11, 12, tzinfo=tzinfo)]
        series = Series(data=values, index=index)

        assert series.index.tz == tzinfo

        # it works! #2443
        repr(series.index[0])

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_tz_aware_asfreq(self, tz):
        dr = date_range('2011-12-01', '2012-07-20', freq='D', tz=tz)

        ser = Series(np.random.randn(len(dr)), index=dr)

        # it works!
        ser.asfreq('T')

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_string_index_alias_tz_aware(self, tz):
        rng = date_range('1/1/2000', periods=10, tz=tz)
        ser = Series(np.random.randn(len(rng)), index=rng)

        result = ser['1/3/2000']
        tm.assert_almost_equal(result, ser[2])

    # TODO: De-duplicate with test below
    def test_series_add_tz_mismatch_converts_to_utc_duplicate(self):
        rng = date_range('1/1/2011', periods=10, freq='H', tz='US/Eastern')
        ser = Series(np.random.randn(len(rng)), index=rng)

        ts_moscow = ser.tz_convert('Europe/Moscow')

        result = ser + ts_moscow
        assert result.index.tz is pytz.utc

        result = ts_moscow + ser
        assert result.index.tz is pytz.utc

    def test_series_add_tz_mismatch_converts_to_utc(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        perm = np.random.permutation(100)[:90]
        ser1 = Series(np.random.randn(90),
                      index=rng.take(perm).tz_convert('US/Eastern'))

        perm = np.random.permutation(100)[:90]
        ser2 = Series(np.random.randn(90),
                      index=rng.take(perm).tz_convert('Europe/Berlin'))

        result = ser1 + ser2

        uts1 = ser1.tz_convert('utc')
        uts2 = ser2.tz_convert('utc')
        expected = uts1 + uts2

        assert result.index.tz == pytz.UTC
        tm.assert_series_equal(result, expected)

    def test_series_add_aware_naive_raises(self):
        rng = date_range('1/1/2011', periods=10, freq='H')
        ser = Series(np.random.randn(len(rng)), index=rng)

        ser_utc = ser.tz_localize('utc')

        with pytest.raises(Exception):
            ser + ser_utc

        with pytest.raises(Exception):
            ser_utc + ser

    def test_series_align_aware(self):
        idx1 = date_range('2001', periods=5, freq='H', tz='US/Eastern')
        ser = Series(np.random.randn(len(idx1)), index=idx1)
        ser_central = ser.tz_convert('US/Central')
        # # different timezones convert to UTC

        new1, new2 = ser.align(ser_central)
        assert new1.index.tz == pytz.UTC
        assert new2.index.tz == pytz.UTC

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_localized_at_time_between_time(self, tzstr):
        from datetime import time
        tz = timezones.maybe_get_tz(tzstr)

        rng = date_range('4/16/2012', '5/1/2012', freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_local = ts.tz_localize(tzstr)

        result = ts_local.at_time(time(10, 0))
        expected = ts.at_time(time(10, 0)).tz_localize(tzstr)
        tm.assert_series_equal(result, expected)
        assert timezones.tz_compare(result.index.tz, tz)

        t1, t2 = time(10, 0), time(11, 0)
        result = ts_local.between_time(t1, t2)
        expected = ts.between_time(t1, t2).tz_localize(tzstr)
        tm.assert_series_equal(result, expected)
        assert timezones.tz_compare(result.index.tz, tz)

    @pytest.mark.parametrize('tzstr', ['Europe/Berlin',
                                       'dateutil/Europe/Berlin'])
    def test_getitem_pydatetime_tz(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        index = date_range(start='2012-12-24 16:00', end='2012-12-24 18:00',
                           freq='H', tz=tzstr)
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp('2012-12-24 17:00', tz=tzstr)

        dt = datetime(2012, 12, 24, 17, 0)
        time_datetime = tslib._localize_pydatetime(dt, tz)
        assert ts[time_pandas] == ts[time_datetime]
