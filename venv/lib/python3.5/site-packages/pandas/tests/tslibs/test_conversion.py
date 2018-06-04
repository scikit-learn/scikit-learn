# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas.util.testing as tm
from pandas import date_range
from pandas._libs.tslib import iNaT
from pandas._libs.tslibs import conversion, timezones


def compare_utc_to_local(tz_didx, utc_didx):
    f = lambda x: conversion.tz_convert_single(x, 'UTC', tz_didx.tz)
    result = conversion.tz_convert(tz_didx.asi8, 'UTC', tz_didx.tz)
    result_single = np.vectorize(f)(tz_didx.asi8)
    tm.assert_numpy_array_equal(result, result_single)


def compare_local_to_utc(tz_didx, utc_didx):
    f = lambda x: conversion.tz_convert_single(x, tz_didx.tz, 'UTC')
    result = conversion.tz_convert(utc_didx.asi8, tz_didx.tz, 'UTC')
    result_single = np.vectorize(f)(utc_didx.asi8)
    tm.assert_numpy_array_equal(result, result_single)


class TestTZConvert(object):

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'Europe/Moscow'])
    def test_tz_convert_single_matches_tz_convert_hourly(self, tz):
        # US: 2014-03-09 - 2014-11-11
        # MOSCOW: 2014-10-26  /  2014-12-31
        tz_didx = date_range('2014-03-01', '2015-01-10', freq='H', tz=tz)
        utc_didx = date_range('2014-03-01', '2015-01-10', freq='H')
        compare_utc_to_local(tz_didx, utc_didx)

        # local tz to UTC can be differ in hourly (or higher) freqs because
        # of DST
        compare_local_to_utc(tz_didx, utc_didx)

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'Europe/Moscow'])
    @pytest.mark.parametrize('freq', ['D', 'A'])
    def test_tz_convert_single_matches_tz_convert(self, tz, freq):
        tz_didx = date_range('2000-01-01', '2020-01-01', freq=freq, tz=tz)
        utc_didx = date_range('2000-01-01', '2020-01-01', freq=freq)
        compare_utc_to_local(tz_didx, utc_didx)
        compare_local_to_utc(tz_didx, utc_didx)

    @pytest.mark.parametrize('arr', [
        pytest.param(np.array([], dtype=np.int64), id='empty'),
        pytest.param(np.array([iNaT], dtype=np.int64), id='all_nat')])
    def test_tz_convert_corner(self, arr):
        result = conversion.tz_convert(arr,
                                       timezones.maybe_get_tz('US/Eastern'),
                                       timezones.maybe_get_tz('Asia/Tokyo'))
        tm.assert_numpy_array_equal(result, arr)
