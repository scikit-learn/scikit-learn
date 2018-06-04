import pytest
import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas.tseries.offsets import Day, Second
from pandas import to_timedelta, timedelta_range


class TestTimedeltas(object):

    def test_timedelta_range(self):

        expected = to_timedelta(np.arange(5), unit='D')
        result = timedelta_range('0 days', periods=5, freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(11), unit='D')
        result = timedelta_range('0 days', '10 days', freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(5), unit='D') + Second(2) + Day()
        result = timedelta_range('1 days, 00:00:02', '5 days, 00:00:02',
                                 freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta([1, 3, 5, 7, 9], unit='D') + Second(2)
        result = timedelta_range('1 days, 00:00:02', periods=5, freq='2D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(50), unit='T') * 30
        result = timedelta_range('0 days', freq='30T', periods=50)
        tm.assert_index_equal(result, expected)

        # GH 11776
        arr = np.arange(10).reshape(2, 5)
        df = pd.DataFrame(np.arange(10).reshape(2, 5))
        for arg in (arr, df):
            with tm.assert_raises_regex(TypeError, "1-d array"):
                to_timedelta(arg)
            for errors in ['ignore', 'raise', 'coerce']:
                with tm.assert_raises_regex(TypeError, "1-d array"):
                    to_timedelta(arg, errors=errors)

        # issue10583
        df = pd.DataFrame(np.random.normal(size=(10, 4)))
        df.index = pd.timedelta_range(start='0s', periods=10, freq='s')
        expected = df.loc[pd.Timedelta('0s'):, :]
        result = df.loc['0s':, :]
        tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize('periods, freq', [
        (3, '2D'), (5, 'D'), (6, '19H12T'), (7, '16H'), (9, '12H')])
    def test_linspace_behavior(self, periods, freq):
        # GH 20976
        result = timedelta_range(start='0 days', end='4 days', periods=periods)
        expected = timedelta_range(start='0 days', end='4 days', freq=freq)
        tm.assert_index_equal(result, expected)

    def test_errors(self):
        # not enough params
        msg = ('Of the four parameters: start, end, periods, and freq, '
               'exactly three must be specified')
        with tm.assert_raises_regex(ValueError, msg):
            timedelta_range(start='0 days')

        with tm.assert_raises_regex(ValueError, msg):
            timedelta_range(end='5 days')

        with tm.assert_raises_regex(ValueError, msg):
            timedelta_range(periods=2)

        with tm.assert_raises_regex(ValueError, msg):
            timedelta_range()

        # too many params
        with tm.assert_raises_regex(ValueError, msg):
            timedelta_range(start='0 days', end='5 days', periods=10, freq='H')
