# -*- coding: utf-8 -*-
"""Tests for PeriodIndex behaving like a vectorized Period scalar"""

from pandas import PeriodIndex, date_range
import pandas.util.testing as tm


class TestPeriodIndexOps(object):
    def test_start_time(self):
        index = PeriodIndex(freq='M', start='2016-01-01', end='2016-05-31')
        expected_index = date_range('2016-01-01', end='2016-05-31', freq='MS')
        tm.assert_index_equal(index.start_time, expected_index)

    def test_end_time(self):
        index = PeriodIndex(freq='M', start='2016-01-01', end='2016-05-31')
        expected_index = date_range('2016-01-01', end='2016-05-31', freq='M')
        tm.assert_index_equal(index.end_time, expected_index)
