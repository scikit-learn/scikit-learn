# -*- coding: utf-8 -*-
"""
Tests for DatetimeIndex timezone-related methods
"""
from datetime import datetime, timedelta, tzinfo
from distutils.version import LooseVersion

import pytest
import pytz
import dateutil
from dateutil.tz import gettz, tzlocal
import numpy as np

import pandas.util.testing as tm
import pandas.util._test_decorators as td

import pandas as pd
from pandas._libs import tslib
from pandas._libs.tslibs import timezones
from pandas.compat import lrange, zip, PY3
from pandas import (DatetimeIndex, date_range, bdate_range,
                    Timestamp, isna, to_datetime, Index)


class FixedOffset(tzinfo):
    """Fixed offset in minutes east from UTC."""

    def __init__(self, offset, name):
        self.__offset = timedelta(minutes=offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return timedelta(0)


fixed_off = FixedOffset(-420, '-07:00')
fixed_off_no_name = FixedOffset(-330, None)


class TestDatetimeIndexTimezones(object):
    # -------------------------------------------------------------
    # DatetimeIndex.tz_convert
    def test_tz_convert_nat(self):
        # GH#5546
        dates = [pd.NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Pacific'))
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Eastern'))
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='UTC'))

        dates = ['2010-12-01 00:00', '2010-12-02 00:00', pd.NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Pacific'))
        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 03:00', '2010-12-02 03:00', pd.NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))

        idx = idx + pd.offsets.Hour(5)
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', pd.NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))
        idx = idx.tz_convert('US/Pacific')
        expected = ['2010-12-01 05:00', '2010-12-02 05:00', pd.NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Pacific'))

        idx = idx + np.timedelta64(3, 'h')
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', pd.NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Pacific'))

        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 11:00', '2010-12-02 11:00', pd.NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))

    @pytest.mark.parametrize('prefix', ['', 'dateutil/'])
    def test_dti_tz_convert_compat_timestamp(self, prefix):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        idx = DatetimeIndex(strdates, tz=prefix + 'US/Eastern')

        conv = idx[0].tz_convert(prefix + 'US/Pacific')
        expected = idx.tz_convert(prefix + 'US/Pacific')[0]

        assert conv == expected

    def test_dti_tz_convert_hour_overflow_dst(self):
        # Regression test for:
        # https://github.com/pandas-dev/pandas/issues/13306

        # sorted case US/Eastern -> UTC
        ts = ['2008-05-12 09:50:00',
              '2008-12-12 09:50:35',
              '2009-05-12 09:50:32']
        tt = DatetimeIndex(ts).tz_localize('US/Eastern')
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # sorted case UTC -> US/Eastern
        ts = ['2008-05-12 13:50:00',
              '2008-12-12 14:50:35',
              '2009-05-12 13:50:32']
        tt = DatetimeIndex(ts).tz_localize('UTC')
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case US/Eastern -> UTC
        ts = ['2008-05-12 09:50:00',
              '2008-12-12 09:50:35',
              '2008-05-12 09:50:32']
        tt = DatetimeIndex(ts).tz_localize('US/Eastern')
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case UTC -> US/Eastern
        ts = ['2008-05-12 13:50:00',
              '2008-12-12 14:50:35',
              '2008-05-12 13:50:32']
        tt = DatetimeIndex(ts).tz_localize('UTC')
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_convert_hour_overflow_dst_timestamps(self, tz):
        # Regression test for GH#13306

        # sorted case US/Eastern -> UTC
        ts = [Timestamp('2008-05-12 09:50:00', tz=tz),
              Timestamp('2008-12-12 09:50:35', tz=tz),
              Timestamp('2009-05-12 09:50:32', tz=tz)]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # sorted case UTC -> US/Eastern
        ts = [Timestamp('2008-05-12 13:50:00', tz='UTC'),
              Timestamp('2008-12-12 14:50:35', tz='UTC'),
              Timestamp('2009-05-12 13:50:32', tz='UTC')]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case US/Eastern -> UTC
        ts = [Timestamp('2008-05-12 09:50:00', tz=tz),
              Timestamp('2008-12-12 09:50:35', tz=tz),
              Timestamp('2008-05-12 09:50:32', tz=tz)]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case UTC -> US/Eastern
        ts = [Timestamp('2008-05-12 13:50:00', tz='UTC'),
              Timestamp('2008-12-12 14:50:35', tz='UTC'),
              Timestamp('2008-05-12 13:50:32', tz='UTC')]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

    @pytest.mark.parametrize('freq, n', [('H', 1), ('T', 60), ('S', 3600)])
    def test_dti_tz_convert_trans_pos_plus_1__bug(self, freq, n):
        # Regression test for tslib.tz_convert(vals, tz1, tz2).
        # See https://github.com/pandas-dev/pandas/issues/4496 for details.
        idx = date_range(datetime(2011, 3, 26, 23),
                         datetime(2011, 3, 27, 1), freq=freq)
        idx = idx.tz_localize('UTC')
        idx = idx.tz_convert('Europe/Moscow')

        expected = np.repeat(np.array([3, 4, 5]), np.array([n, n, 1]))
        tm.assert_index_equal(idx.hour, Index(expected))

    def test_dti_tz_convert_dst(self):
        for freq, n in [('H', 1), ('T', 60), ('S', 3600)]:
            # Start DST
            idx = date_range('2014-03-08 23:00', '2014-03-09 09:00', freq=freq,
                             tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([18, 19, 20, 21, 22, 23,
                                           0, 1, 3, 4, 5]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            idx = date_range('2014-03-08 18:00', '2014-03-09 05:00', freq=freq,
                             tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            # End DST
            idx = date_range('2014-11-01 23:00', '2014-11-02 09:00', freq=freq,
                             tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([19, 20, 21, 22, 23,
                                           0, 1, 1, 2, 3, 4]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            idx = date_range('2014-11-01 18:00', '2014-11-02 05:00', freq=freq,
                             tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([22, 23, 0, 1, 2, 3, 4, 5, 6,
                                           7, 8, 9, 10]),
                                 np.array([n, n, n, n, n, n, n, n, n,
                                           n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

        # daily
        # Start DST
        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D',
                         tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx.hour, Index([19, 19]))

        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D',
                         tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx.hour, Index([5, 5]))

        # End DST
        idx = date_range('2014-11-01 00:00', '2014-11-02 00:00', freq='D',
                         tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx.hour, Index([20, 20]))

        idx = date_range('2014-11-01 00:00', '2014-11-02 000:00', freq='D',
                         tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx.hour, Index([4, 4]))

    def test_tz_convert_roundtrip(self, tz_aware_fixture):
        tz = tz_aware_fixture
        idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M',
                          tz='UTC')
        exp1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')

        idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D',
                          tz='UTC')
        exp2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')

        idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H',
                          tz='UTC')
        exp3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')

        idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T',
                          tz='UTC')
        exp4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')

        for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3),
                              (idx4, exp4)]:
            converted = idx.tz_convert(tz)
            reset = converted.tz_convert(None)
            tm.assert_index_equal(reset, expected)
            assert reset.tzinfo is None
            expected = converted.tz_convert('UTC').tz_localize(None)
            tm.assert_index_equal(reset, expected)

    def test_dti_tz_convert_tzlocal(self):
        # GH#13583
        # tz_convert doesn't affect to internal
        dti = date_range(start='2001-01-01', end='2001-03-01', tz='UTC')
        dti2 = dti.tz_convert(dateutil.tz.tzlocal())
        tm.assert_numpy_array_equal(dti2.asi8, dti.asi8)

        dti = date_range(start='2001-01-01', end='2001-03-01',
                         tz=dateutil.tz.tzlocal())
        dti2 = dti.tz_convert(None)
        tm.assert_numpy_array_equal(dti2.asi8, dti.asi8)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern',
                                    pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_convert_utc_to_local_no_modify(self, tz):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(tz)

        # Values are unmodified
        tm.assert_numpy_array_equal(rng.asi8, rng_eastern.asi8)

        assert timezones.tz_compare(rng_eastern.tz, timezones.maybe_get_tz(tz))

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_tz_convert_unsorted(self, tzstr):
        dr = date_range('2012-03-09', freq='H', periods=100, tz='utc')
        dr = dr.tz_convert(tzstr)

        result = dr[::-1].hour
        exp = dr.hour[::-1]
        tm.assert_almost_equal(result, exp)

    # -------------------------------------------------------------
    # DatetimeIndex.tz_localize

    def test_dti_tz_localize_nonexistent_raise_coerce(self):
        # GH#13057
        times = ['2015-03-08 01:00', '2015-03-08 02:00', '2015-03-08 03:00']
        index = DatetimeIndex(times)
        tz = 'US/Eastern'
        with pytest.raises(pytz.NonExistentTimeError):
            index.tz_localize(tz=tz)

        with pytest.raises(pytz.NonExistentTimeError):
            index.tz_localize(tz=tz, errors='raise')

        result = index.tz_localize(tz=tz, errors='coerce')
        test_times = ['2015-03-08 01:00-05:00', 'NaT',
                      '2015-03-08 03:00-04:00']
        dti = DatetimeIndex(test_times)
        expected = dti.tz_localize('UTC').tz_convert('US/Eastern')
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_infer(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour
        # With no repeated hours, we cannot infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.AmbiguousTimeError):
            dr.tz_localize(tz)

        # With repeated hours, we can infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']
        di = DatetimeIndex(times)
        localized = di.tz_localize(tz, ambiguous='infer')
        tm.assert_index_equal(dr, localized)
        tm.assert_index_equal(dr, DatetimeIndex(times, tz=tz,
                                                ambiguous='infer'))

        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=pd.offsets.Hour())
        localized = dr.tz_localize(tz)
        localized_infer = dr.tz_localize(tz, ambiguous='infer')
        tm.assert_index_equal(localized, localized_infer)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_times(self, tz):
        # March 13, 2011, spring forward, skip from 2 AM to 3 AM
        dr = date_range(datetime(2011, 3, 13, 1, 30), periods=3,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.NonExistentTimeError):
            dr.tz_localize(tz)

        # after dst transition, it works
        dr = date_range(datetime(2011, 3, 13, 3, 30), periods=3,
                        freq=pd.offsets.Hour(), tz=tz)

        # November 6, 2011, fall back, repeat 2 AM hour
        dr = date_range(datetime(2011, 11, 6, 1, 30), periods=3,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.AmbiguousTimeError):
            dr.tz_localize(tz)

        # UTC is OK
        dr = date_range(datetime(2011, 3, 13), periods=48,
                        freq=pd.offsets.Minute(30), tz=pytz.utc)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_pass_dates_to_utc(self, tzstr):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']

        idx = DatetimeIndex(strdates)
        conv = idx.tz_localize(tzstr)

        fromdates = DatetimeIndex(strdates, tz=tzstr)

        assert conv.tz == fromdates.tz
        tm.assert_numpy_array_equal(conv.values, fromdates.values)

    @pytest.mark.parametrize('prefix', ['', 'dateutil/'])
    def test_dti_tz_localize(self, prefix):
        tzstr = prefix + 'US/Eastern'
        dti = DatetimeIndex(start='1/1/2005', end='1/1/2005 0:00:30.256',
                            freq='L')
        dti2 = dti.tz_localize(tzstr)

        dti_utc = DatetimeIndex(start='1/1/2005 05:00',
                                end='1/1/2005 5:00:30.256', freq='L', tz='utc')

        tm.assert_numpy_array_equal(dti2.values, dti_utc.values)

        dti3 = dti2.tz_convert(prefix + 'US/Pacific')
        tm.assert_numpy_array_equal(dti3.values, dti_utc.values)

        dti = DatetimeIndex(start='11/6/2011 1:59', end='11/6/2011 2:00',
                            freq='L')
        with pytest.raises(pytz.AmbiguousTimeError):
            dti.tz_localize(tzstr)

        dti = DatetimeIndex(start='3/13/2011 1:59', end='3/13/2011 2:00',
                            freq='L')
        with pytest.raises(pytz.NonExistentTimeError):
            dti.tz_localize(tzstr)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern',
                                    pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_utc_conversion(self, tz):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range('3/10/2012', '3/11/2012', freq='30T')

        converted = rng.tz_localize(tz)
        expected_naive = rng + pd.offsets.Hour(5)
        tm.assert_numpy_array_equal(converted.asi8, expected_naive.asi8)

        # DST ambiguity, this should fail
        rng = date_range('3/11/2012', '3/12/2012', freq='30T')
        # Is this really how it should fail??
        with pytest.raises(pytz.NonExistentTimeError):
            rng.tz_localize(tz)

    def test_dti_tz_localize_roundtrip(self, tz_aware_fixture):
        tz = tz_aware_fixture

        idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')
        idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')
        idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')
        idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')
        for idx in [idx1, idx2, idx3, idx4]:
            localized = idx.tz_localize(tz)
            expected = date_range(start=idx[0], end=idx[-1], freq=idx.freq,
                                  tz=tz)
            tm.assert_index_equal(localized, expected)
            with pytest.raises(TypeError):
                localized.tz_localize(tz)

            reset = localized.tz_localize(None)
            tm.assert_index_equal(reset, idx)
            assert reset.tzinfo is None

    def test_dti_tz_localize_naive(self):
        rng = date_range('1/1/2011', periods=100, freq='H')

        conv = rng.tz_localize('US/Pacific')
        exp = date_range('1/1/2011', periods=100, freq='H', tz='US/Pacific')

        tm.assert_index_equal(conv, exp)

    def test_dti_tz_localize_tzlocal(self):
        # GH#13583
        offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
        offset = int(offset.total_seconds() * 1000000000)

        dti = date_range(start='2001-01-01', end='2001-03-01')
        dti2 = dti.tz_localize(dateutil.tz.tzlocal())
        tm.assert_numpy_array_equal(dti2.asi8 + offset, dti.asi8)

        dti = date_range(start='2001-01-01', end='2001-03-01',
                         tz=dateutil.tz.tzlocal())
        dti2 = dti.tz_localize(None)
        tm.assert_numpy_array_equal(dti2.asi8 - offset, dti.asi8)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_nat(self, tz):
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']
        di = DatetimeIndex(times)
        localized = di.tz_localize(tz, ambiguous='NaT')

        times = ['11/06/2011 00:00', np.NaN, np.NaN, '11/06/2011 02:00',
                 '11/06/2011 03:00']
        di_test = DatetimeIndex(times, tz='US/Eastern')

        # left dtype is datetime64[ns, US/Eastern]
        # right is datetime64[ns, tzfile('/usr/share/zoneinfo/US/Eastern')]
        tm.assert_numpy_array_equal(di_test.values, localized.values)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_flags(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour

        # Pass in flags to determine right dst transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']

        # Test tz_localize
        di = DatetimeIndex(times)
        is_dst = [1, 1, 0, 0, 0]
        localized = di.tz_localize(tz, ambiguous=is_dst)
        tm.assert_index_equal(dr, localized)
        tm.assert_index_equal(dr, DatetimeIndex(times, tz=tz,
                                                ambiguous=is_dst))

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst))
        tm.assert_index_equal(dr, localized)

        localized = di.tz_localize(tz,
                                   ambiguous=np.array(is_dst).astype('bool'))
        tm.assert_index_equal(dr, localized)

        # Test constructor
        localized = DatetimeIndex(times, tz=tz, ambiguous=is_dst)
        tm.assert_index_equal(dr, localized)

        # Test duplicate times where inferring the dst fails
        times += times
        di = DatetimeIndex(times)

        # When the sizes are incompatible, make sure error is raised
        with pytest.raises(Exception):
            di.tz_localize(tz, ambiguous=is_dst)

        # When sizes are compatible and there are repeats ('infer' won't work)
        is_dst = np.hstack((is_dst, is_dst))
        localized = di.tz_localize(tz, ambiguous=is_dst)
        dr = dr.append(dr)
        tm.assert_index_equal(dr, localized)

        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=pd.offsets.Hour())
        is_dst = np.array([1] * 10)
        localized = dr.tz_localize(tz)
        localized_is_dst = dr.tz_localize(tz, ambiguous=is_dst)
        tm.assert_index_equal(localized, localized_is_dst)

    # TODO: belongs outside tz_localize tests?
    @pytest.mark.parametrize('tz', ['Europe/London', 'dateutil/Europe/London'])
    def test_dti_construction_ambiguous_endpoint(self, tz):
        # construction with an ambiguous end-point
        # GH#11626

        # FIXME: This next block fails to raise; it was taken from an older
        # version of this test that had an indention mistake that caused it
        # to not get executed.
        # with pytest.raises(pytz.AmbiguousTimeError):
        #    date_range("2013-10-26 23:00", "2013-10-27 01:00",
        #               tz="Europe/London", freq="H")

        times = date_range("2013-10-26 23:00", "2013-10-27 01:00", freq="H",
                           tz=tz, ambiguous='infer')
        assert times[0] == Timestamp('2013-10-26 23:00', tz=tz, freq="H")

        if str(tz).startswith('dateutil'):
            if LooseVersion(dateutil.__version__) < LooseVersion('2.6.0'):
                # see GH#14621
                assert times[-1] == Timestamp('2013-10-27 01:00:00+0000',
                                              tz=tz, freq="H")
            elif LooseVersion(dateutil.__version__) > LooseVersion('2.6.0'):
                # fixed ambiguous behavior
                assert times[-1] == Timestamp('2013-10-27 01:00:00+0100',
                                              tz=tz, freq="H")
        else:
            assert times[-1] == Timestamp('2013-10-27 01:00:00+0000',
                                          tz=tz, freq="H")

    def test_dti_tz_localize_bdate_range(self):
        dr = pd.bdate_range('1/1/2009', '1/1/2010')
        dr_utc = pd.bdate_range('1/1/2009', '1/1/2010', tz=pytz.utc)
        localized = dr.tz_localize(pytz.utc)
        tm.assert_index_equal(dr_utc, localized)

    # -------------------------------------------------------------
    # DatetimeIndex.normalize

    def test_normalize_tz(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='US/Eastern')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

        rng = date_range('1/1/2000 9:30', periods=10, freq='D', tz='UTC')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D', tz='UTC')
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

        rng = date_range('1/1/2000 9:30', periods=10, freq='D', tz=tzlocal())
        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D', tz=tzlocal())
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

    @td.skip_if_windows
    @pytest.mark.parametrize('timezone', ['US/Pacific', 'US/Eastern', 'UTC',
                                          'Asia/Kolkata', 'Asia/Shanghai',
                                          'Australia/Canberra'])
    def test_normalize_tz_local(self, timezone):
        # GH#13459
        with tm.set_timezone(timezone):
            rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                             tz=tzlocal())

            result = rng.normalize()
            expected = date_range('1/1/2000', periods=10, freq='D',
                                  tz=tzlocal())
            tm.assert_index_equal(result, expected)

            assert result.is_normalized
            assert not rng.is_normalized

    # ------------------------------------------------------------
    # DatetimeIndex.__new__

    @pytest.mark.parametrize('prefix', ['', 'dateutil/'])
    def test_dti_constructor_static_tzinfo(self, prefix):
        # it works!
        index = DatetimeIndex([datetime(2012, 1, 1)], tz=prefix + 'EST')
        index.hour
        index[0]

    def test_dti_constructor_with_fixed_tz(self):
        off = FixedOffset(420, '+07:00')
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        assert off == rng.tz

        rng2 = date_range(start, periods=len(rng), tz=off)
        tm.assert_index_equal(rng, rng2)

        rng3 = date_range('3/11/2012 05:00:00+07:00',
                          '6/11/2012 05:00:00+07:00')
        assert (rng.values == rng3.values).all()

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_convert_datetime_list(self, tzstr):
        dr = date_range('2012-06-02', periods=10,
                        tz=tzstr, name='foo')
        dr2 = DatetimeIndex(list(dr), name='foo')
        tm.assert_index_equal(dr, dr2)
        assert dr.tz == dr2.tz
        assert dr2.name == 'foo'

    def test_dti_construction_univalent(self):
        rng = date_range('03/12/2012 00:00', periods=10, freq='W-FRI',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(data=rng, tz='US/Eastern')
        tm.assert_index_equal(rng, rng2)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_from_tzaware_datetime(self, tz):
        d = [datetime(2012, 8, 19, tzinfo=tz)]

        index = DatetimeIndex(d)
        assert timezones.tz_compare(index.tz, tz)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_constructors(self, tzstr):
        """ Test different DatetimeIndex constructions with timezone
        Follow-up of GH#4229
        """

        arr = ['11/10/2005 08:00:00', '11/10/2005 09:00:00']

        idx1 = to_datetime(arr).tz_localize(tzstr)
        idx2 = DatetimeIndex(start="2005-11-10 08:00:00", freq='H', periods=2,
                             tz=tzstr)
        idx3 = DatetimeIndex(arr, tz=tzstr)
        idx4 = DatetimeIndex(np.array(arr), tz=tzstr)

        for other in [idx2, idx3, idx4]:
            tm.assert_index_equal(idx1, other)

    # -------------------------------------------------------------
    # Unsorted

    def test_join_utc_convert(self, join_type):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng.tz_convert('US/Eastern')
        right = rng.tz_convert('Europe/Berlin')

        result = left.join(left[:-5], how=join_type)
        assert isinstance(result, DatetimeIndex)
        assert result.tz == left.tz

        result = left.join(right[:-5], how=join_type)
        assert isinstance(result, DatetimeIndex)
        assert result.tz.zone == 'UTC'

    def test_dti_drop_dont_lose_tz(self):
        # GH#2621
        ind = date_range("2012-12-01", periods=10, tz="utc")
        ind = ind.drop(ind[-1])

        assert ind.tz is not None

    def test_date_range_localize(self):
        rng = date_range('3/11/2012 03:00', periods=15, freq='H',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 03:00', '3/11/2012 04:00'],
                             tz='US/Eastern')
        rng3 = date_range('3/11/2012 03:00', periods=15, freq='H')
        rng3 = rng3.tz_localize('US/Eastern')

        tm.assert_index_equal(rng, rng3)

        # DST transition time
        val = rng[0]
        exp = Timestamp('3/11/2012 03:00', tz='US/Eastern')

        assert val.hour == 3
        assert exp.hour == 3
        assert val == exp  # same UTC value
        tm.assert_index_equal(rng[:2], rng2)

        # Right before the DST transition
        rng = date_range('3/11/2012 00:00', periods=2, freq='H',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 00:00', '3/11/2012 01:00'],
                             tz='US/Eastern')
        tm.assert_index_equal(rng, rng2)
        exp = Timestamp('3/11/2012 00:00', tz='US/Eastern')
        assert exp.hour == 0
        assert rng[0] == exp
        exp = Timestamp('3/11/2012 01:00', tz='US/Eastern')
        assert exp.hour == 1
        assert rng[1] == exp

        rng = date_range('3/11/2012 00:00', periods=10, freq='H',
                         tz='US/Eastern')
        assert rng[2].hour == 3

    def test_timestamp_equality_different_timezones(self):
        utc_range = date_range('1/1/2000', periods=20, tz='UTC')
        eastern_range = utc_range.tz_convert('US/Eastern')
        berlin_range = utc_range.tz_convert('Europe/Berlin')

        for a, b, c in zip(utc_range, eastern_range, berlin_range):
            assert a == b
            assert b == c
            assert a == c

        assert (utc_range == eastern_range).all()
        assert (utc_range == berlin_range).all()
        assert (berlin_range == eastern_range).all()

    def test_dti_intersection(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng[10:90][::-1]
        right = rng[20:80][::-1]

        assert left.tz == rng.tz
        result = left.intersection(right)
        assert result.tz == left.tz

    def test_dti_equals_with_tz(self):
        left = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        right = date_range('1/1/2011', periods=100, freq='H', tz='US/Eastern')

        assert not left.equals(right)

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_nat(self, tzstr):
        idx = DatetimeIndex([Timestamp("2013-1-1", tz=tzstr), pd.NaT])

        assert isna(idx[1])
        assert idx[0].tzinfo is not None

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_astype_asobject_tzinfos(self, tzstr):
        # GH#1345

        # dates around a dst transition
        rng = date_range('2/13/2010', '5/6/2010', tz=tzstr)

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            assert x == exval
            assert x.tzinfo == exval.tzinfo

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            assert x == exval
            assert x.tzinfo == exval.tzinfo

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_with_timezone_repr(self, tzstr):
        rng = date_range('4/13/2010', '5/6/2010')

        rng_eastern = rng.tz_localize(tzstr)

        rng_repr = repr(rng_eastern)
        assert '2010-04-13 00:00:00' in rng_repr

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_take_dont_lose_meta(self, tzstr):
        rng = date_range('1/1/2000', periods=20, tz=tzstr)

        result = rng.take(lrange(5))
        assert result.tz == rng.tz
        assert result.freq == rng.freq

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_utc_box_timestamp_and_localize(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(tzstr)

        expected = rng[-1].astimezone(tz)

        stamp = rng_eastern[-1]
        assert stamp == expected
        assert stamp.tzinfo == expected.tzinfo

        # right tzinfo
        rng = date_range('3/13/2012', '3/14/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(tzstr)
        # test not valid for dateutil timezones.
        # assert 'EDT' in repr(rng_eastern[0].tzinfo)
        assert ('EDT' in repr(rng_eastern[0].tzinfo) or
                'tzfile' in repr(rng_eastern[0].tzinfo))

    def test_dti_to_pydatetime(self):
        dt = dateutil.parser.parse('2012-06-13T01:39:00Z')
        dt = dt.replace(tzinfo=tzlocal())

        arr = np.array([dt], dtype=object)

        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc

        rng = date_range('2012-11-03 03:00', '2012-11-05 03:00', tz=tzlocal())
        arr = rng.to_pydatetime()
        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc

    def test_dti_to_pydatetime_fizedtz(self):
        dates = np.array([datetime(2000, 1, 1, tzinfo=fixed_off),
                          datetime(2000, 1, 2, tzinfo=fixed_off),
                          datetime(2000, 1, 3, tzinfo=fixed_off)])
        dti = DatetimeIndex(dates)

        result = dti.to_pydatetime()
        tm.assert_numpy_array_equal(dates, result)

        result = dti._mpl_repr()
        tm.assert_numpy_array_equal(dates, result)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Central'),
                                    gettz('US/Central')])
    def test_with_tz(self, tz):
        # just want it to work
        start = datetime(2011, 3, 12, tzinfo=pytz.utc)
        dr = bdate_range(start, periods=50, freq=pd.offsets.Hour())
        assert dr.tz is pytz.utc

        # DateRange with naive datetimes
        dr = bdate_range('1/1/2005', '1/1/2009', tz=pytz.utc)
        dr = bdate_range('1/1/2005', '1/1/2009', tz=tz)

        # normalized
        central = dr.tz_convert(tz)
        assert central.tz is tz
        naive = central[0].to_pydatetime().replace(tzinfo=None)
        comp = tslib._localize_pydatetime(naive, tz).tzinfo
        assert central[0].tz is comp

        # compare vs a localized tz
        naive = dr[0].to_pydatetime().replace(tzinfo=None)
        comp = tslib._localize_pydatetime(naive, tz).tzinfo
        assert central[0].tz is comp

        # datetimes with tzinfo set
        dr = bdate_range(datetime(2005, 1, 1, tzinfo=pytz.utc),
                         datetime(2009, 1, 1, tzinfo=pytz.utc))
        with pytest.raises(Exception):
            bdate_range(datetime(2005, 1, 1, tzinfo=pytz.utc), '1/1/2009',
                        tz=tz)

    @pytest.mark.parametrize('prefix', ['', 'dateutil/'])
    def test_field_access_localize(self, prefix):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        rng = DatetimeIndex(strdates, tz=prefix + 'US/Eastern')
        assert (rng.hour == 0).all()

        # a more unusual time zone, #1946
        dr = date_range('2011-10-02 00:00', freq='h', periods=10,
                        tz=prefix + 'America/Atikokan')

        expected = Index(np.arange(10, dtype=np.int64))
        tm.assert_index_equal(dr.hour, expected)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_convert_tz_aware_datetime_datetime(self, tz):
        # GH#1581
        dates = [datetime(2000, 1, 1), datetime(2000, 1, 2),
                 datetime(2000, 1, 3)]

        dates_aware = [tslib._localize_pydatetime(x, tz) for x in dates]
        result = DatetimeIndex(dates_aware)
        assert timezones.tz_compare(result.tz, tz)

        converted = to_datetime(dates_aware, utc=True)
        ex_vals = np.array([Timestamp(x).value for x in dates_aware])
        tm.assert_numpy_array_equal(converted.asi8, ex_vals)
        assert converted.tz is pytz.utc

    def test_dti_union_aware(self):
        # non-overlapping
        rng = date_range("2012-11-15 00:00:00", periods=6, freq="H",
                         tz="US/Central")

        rng2 = date_range("2012-11-15 12:00:00", periods=6, freq="H",
                          tz="US/Eastern")

        result = rng.union(rng2)
        assert result.tz.zone == 'UTC'

    @pytest.mark.parametrize('tz', [None, 'UTC', "US/Central",
                                    dateutil.tz.tzoffset(None, -28800)])
    @pytest.mark.usefixtures("datetime_tz_utc")
    @pytest.mark.skipif(not PY3, reason="datetime.timezone not in PY2")
    def test_iteration_preserves_nanoseconds(self, tz):
        # GH 19603
        index = DatetimeIndex(["2018-02-08 15:00:00.168456358",
                               "2018-02-08 15:00:00.168456359"], tz=tz)
        for i, ts in enumerate(index):
            assert ts == index[i]


class TestDateRange(object):
    """Tests for date_range with timezones"""
    def test_hongkong_tz_convert(self):
        # GH#1673 smoke test
        dr = date_range('2012-01-01', '2012-01-10', freq='D', tz='Hongkong')

        # it works!
        dr.hour

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_date_range_span_dst_transition(self, tzstr):
        # GH#1778

        # Standard -> Daylight Savings Time
        dr = date_range('03/06/2012 00:00', periods=200, freq='W-FRI',
                        tz='US/Eastern')

        assert (dr.hour == 0).all()

        dr = date_range('2012-11-02', periods=10, tz=tzstr)
        assert (dr.hour == 0).all()

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_date_range_timezone_str_argument(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)
        result = date_range('1/1/2000', periods=10, tz=tzstr)
        expected = date_range('1/1/2000', periods=10, tz=tz)

        tm.assert_index_equal(result, expected)

    def test_date_range_with_fixedoffset_noname(self):
        off = fixed_off_no_name
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        assert off == rng.tz

        idx = Index([start, end])
        assert off == idx.tz

    @pytest.mark.parametrize('tzstr', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_date_range_with_tz(self, tzstr):
        stamp = Timestamp('3/11/2012 05:00', tz=tzstr)
        assert stamp.hour == 5

        rng = date_range('3/11/2012 04:00', periods=10, freq='H',
                         tz=tzstr)

        assert stamp == rng[1]


class TestToDatetime(object):
    """Tests for the to_datetime constructor with timezones"""
    def test_to_datetime_utc(self):
        arr = np.array([dateutil.parser.parse('2012-06-13T01:39:00Z')],
                       dtype=object)

        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc

    def test_to_datetime_fixed_offset(self):
        dates = [datetime(2000, 1, 1, tzinfo=fixed_off),
                 datetime(2000, 1, 2, tzinfo=fixed_off),
                 datetime(2000, 1, 3, tzinfo=fixed_off)]
        result = to_datetime(dates)
        assert result.tz == fixed_off
