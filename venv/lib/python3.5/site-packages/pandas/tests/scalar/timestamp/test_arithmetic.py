# -*- coding: utf-8 -*-
from datetime import datetime, timedelta

import pytest
import numpy as np

from pandas.compat import long
from pandas.tseries import offsets
from pandas import Timestamp, Timedelta


class TestTimestampArithmetic(object):
    def test_overflow_offset(self):
        # xref https://github.com/statsmodels/statsmodels/issues/3374
        # ends up multiplying really large numbers which overflow

        stamp = Timestamp('2017-01-13 00:00:00', freq='D')
        offset = 20169940 * offsets.Day(1)

        with pytest.raises(OverflowError):
            stamp + offset

        with pytest.raises(OverflowError):
            offset + stamp

        with pytest.raises(OverflowError):
            stamp - offset

    def test_delta_preserve_nanos(self):
        val = Timestamp(long(1337299200000000123))
        result = val + timedelta(1)
        assert result.nanosecond == val.nanosecond

    def test_timestamp_sub_datetime(self):
        dt = datetime(2013, 10, 12)
        ts = Timestamp(datetime(2013, 10, 13))
        assert (ts - dt).days == 1
        assert (dt - ts).days == -1

    def test_addition_subtraction_types(self):
        # Assert on the types resulting from Timestamp +/- various date/time
        # objects
        dt = datetime(2014, 3, 4)
        td = timedelta(seconds=1)
        # build a timestamp with a frequency, since then it supports
        # addition/subtraction of integers
        ts = Timestamp(dt, freq='D')

        assert type(ts + 1) == Timestamp
        assert type(ts - 1) == Timestamp

        # Timestamp + datetime not supported, though subtraction is supported
        # and yields timedelta more tests in tseries/base/tests/test_base.py
        assert type(ts - dt) == Timedelta
        assert type(ts + td) == Timestamp
        assert type(ts - td) == Timestamp

        # Timestamp +/- datetime64 not supported, so not tested (could possibly
        # assert error raised?)
        td64 = np.timedelta64(1, 'D')
        assert type(ts + td64) == Timestamp
        assert type(ts - td64) == Timestamp

    def test_addition_subtraction_preserve_frequency(self):
        ts = Timestamp('2014-03-05', freq='D')
        td = timedelta(days=1)
        original_freq = ts.freq

        assert (ts + 1).freq == original_freq
        assert (ts - 1).freq == original_freq
        assert (ts + td).freq == original_freq
        assert (ts - td).freq == original_freq

        td64 = np.timedelta64(1, 'D')
        assert (ts + td64).freq == original_freq
        assert (ts - td64).freq == original_freq
