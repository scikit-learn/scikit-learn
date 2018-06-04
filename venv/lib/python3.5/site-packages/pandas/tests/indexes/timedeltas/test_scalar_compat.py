# -*- coding: utf-8 -*-
"""
Tests for TimedeltaIndex methods behaving like their Timedelta counterparts
"""

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import timedelta_range, Timedelta, TimedeltaIndex, Index, Series


class TestVectorizedTimedelta(object):
    def test_tdi_total_seconds(self):
        # GH#10939
        # test index
        rng = timedelta_range('1 days, 10:11:12.100123456', periods=2,
                              freq='s')
        expt = [1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9,
                1 * 86400 + 10 * 3600 + 11 * 60 + 13 + 100123456. / 1e9]
        tm.assert_almost_equal(rng.total_seconds(), Index(expt))

        # test Series
        ser = Series(rng)
        s_expt = Series(expt, index=[0, 1])
        tm.assert_series_equal(ser.dt.total_seconds(), s_expt)

        # with nat
        ser[1] = np.nan
        s_expt = Series([1 * 86400 + 10 * 3600 + 11 * 60 +
                         12 + 100123456. / 1e9, np.nan], index=[0, 1])
        tm.assert_series_equal(ser.dt.total_seconds(), s_expt)

        # with both nat
        ser = Series([np.nan, np.nan], dtype='timedelta64[ns]')
        tm.assert_series_equal(ser.dt.total_seconds(),
                               Series([np.nan, np.nan], index=[0, 1]))

    def test_tdi_round(self):
        td = pd.timedelta_range(start='16801 days', periods=5, freq='30Min')
        elt = td[1]

        expected_rng = TimedeltaIndex([Timedelta('16801 days 00:00:00'),
                                       Timedelta('16801 days 00:00:00'),
                                       Timedelta('16801 days 01:00:00'),
                                       Timedelta('16801 days 02:00:00'),
                                       Timedelta('16801 days 02:00:00')])
        expected_elt = expected_rng[1]

        tm.assert_index_equal(td.round(freq='H'), expected_rng)
        assert elt.round(freq='H') == expected_elt

        msg = pd._libs.tslibs.frequencies._INVALID_FREQ_ERROR
        with tm.assert_raises_regex(ValueError, msg):
            td.round(freq='foo')
        with tm.assert_raises_regex(ValueError, msg):
            elt.round(freq='foo')

        msg = "<MonthEnd> is a non-fixed frequency"
        with tm.assert_raises_regex(ValueError, msg):
            td.round(freq='M')
        with tm.assert_raises_regex(ValueError, msg):
            elt.round(freq='M')
