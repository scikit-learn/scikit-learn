from datetime import timedelta

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import TimedeltaIndex, timedelta_range, compat, Index, Timedelta


class TestGetItem(object):
    def test_getitem(self):
        idx1 = timedelta_range('1 day', '31 day', freq='D', name='idx')

        for idx in [idx1]:
            result = idx[0]
            assert result == Timedelta('1 day')

            result = idx[0:5]
            expected = timedelta_range('1 day', '5 day', freq='D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[0:10:2]
            expected = timedelta_range('1 day', '9 day', freq='2D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[-20:-5:3]
            expected = timedelta_range('12 day', '24 day', freq='3D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[4::-1]
            expected = TimedeltaIndex(['5 day', '4 day', '3 day',
                                       '2 day', '1 day'],
                                      freq='-1D', name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq


class TestWhere(object):
    # placeholder for symmetry with DatetimeIndex and PeriodIndex tests
    pass


class TestTake(object):
    def test_take(self):
        # GH 10295
        idx1 = timedelta_range('1 day', '31 day', freq='D', name='idx')

        for idx in [idx1]:
            result = idx.take([0])
            assert result == Timedelta('1 day')

            result = idx.take([-1])
            assert result == Timedelta('31 day')

            result = idx.take([0, 1, 2])
            expected = timedelta_range('1 day', '3 day', freq='D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([0, 2, 4])
            expected = timedelta_range('1 day', '5 day', freq='2D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([7, 4, 1])
            expected = timedelta_range('8 day', '2 day', freq='-3D',
                                       name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([3, 2, 5])
            expected = TimedeltaIndex(['4 day', '3 day', '6 day'], name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq is None

            result = idx.take([-3, 2, 5])
            expected = TimedeltaIndex(['29 day', '3 day', '6 day'], name='idx')
            tm.assert_index_equal(result, expected)
            assert result.freq is None

    def test_take_invalid_kwargs(self):
        idx = timedelta_range('1 day', '31 day', freq='D', name='idx')
        indices = [1, 6, 5, 9, 10, 13, 15, 3]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        tm.assert_raises_regex(TypeError, msg, idx.take,
                               indices, foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, idx.take,
                               indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, idx.take,
                               indices, mode='clip')

    # TODO: This method came from test_timedelta; de-dup with version above
    def test_take2(self):
        tds = ['1day 02:00:00', '1 day 04:00:00', '1 day 10:00:00']
        idx = TimedeltaIndex(start='1d', end='2d', freq='H', name='idx')
        expected = TimedeltaIndex(tds, freq=None, name='idx')

        taken1 = idx.take([2, 4, 10])
        taken2 = idx[[2, 4, 10]]

        for taken in [taken1, taken2]:
            tm.assert_index_equal(taken, expected)
            assert isinstance(taken, TimedeltaIndex)
            assert taken.freq is None
            assert taken.name == expected.name

    def test_take_fill_value(self):
        # GH 12631
        idx = TimedeltaIndex(['1 days', '2 days', '3 days'],
                             name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = TimedeltaIndex(['2 days', '1 days', '3 days'],
                                  name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = TimedeltaIndex(['2 days', '1 days', 'NaT'],
                                  name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = TimedeltaIndex(['2 days', '1 days', '3 days'],
                                  name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with pytest.raises(IndexError):
            idx.take(np.array([1, -5]))


class TestTimedeltaIndex(object):

    def test_insert(self):

        idx = TimedeltaIndex(['4day', '1day', '2day'], name='idx')

        result = idx.insert(2, timedelta(days=5))
        exp = TimedeltaIndex(['4day', '1day', '5day', '2day'], name='idx')
        tm.assert_index_equal(result, exp)

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, 'inserted')
        expected = Index([Timedelta('4day'), 'inserted', Timedelta('1day'),
                          Timedelta('2day')], name='idx')
        assert not isinstance(result, TimedeltaIndex)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

        idx = timedelta_range('1day 00:00:01', periods=3, freq='s', name='idx')

        # preserve freq
        expected_0 = TimedeltaIndex(['1day', '1day 00:00:01', '1day 00:00:02',
                                     '1day 00:00:03'],
                                    name='idx', freq='s')
        expected_3 = TimedeltaIndex(['1day 00:00:01', '1day 00:00:02',
                                     '1day 00:00:03', '1day 00:00:04'],
                                    name='idx', freq='s')

        # reset freq to None
        expected_1_nofreq = TimedeltaIndex(['1day 00:00:01', '1day 00:00:01',
                                            '1day 00:00:02', '1day 00:00:03'],
                                           name='idx', freq=None)
        expected_3_nofreq = TimedeltaIndex(['1day 00:00:01', '1day 00:00:02',
                                            '1day 00:00:03', '1day 00:00:05'],
                                           name='idx', freq=None)

        cases = [(0, Timedelta('1day'), expected_0),
                 (-3, Timedelta('1day'), expected_0),
                 (3, Timedelta('1day 00:00:04'), expected_3),
                 (1, Timedelta('1day 00:00:01'), expected_1_nofreq),
                 (3, Timedelta('1day 00:00:05'), expected_3_nofreq)]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

        # GH 18295 (test missing)
        expected = TimedeltaIndex(['1day', pd.NaT, '2day', '3day'])
        for na in (np.nan, pd.NaT, None):
            result = timedelta_range('1day', '3day').insert(1, na)
            tm.assert_index_equal(result, expected)

    def test_delete(self):
        idx = timedelta_range(start='1 Days', periods=5, freq='D', name='idx')

        # prserve freq
        expected_0 = timedelta_range(start='2 Days', periods=4, freq='D',
                                     name='idx')
        expected_4 = timedelta_range(start='1 Days', periods=4, freq='D',
                                     name='idx')

        # reset freq to None
        expected_1 = TimedeltaIndex(
            ['1 day', '3 day', '4 day', '5 day'], freq=None, name='idx')

        cases = {0: expected_0,
                 -5: expected_0,
                 -1: expected_4,
                 4: expected_4,
                 1: expected_1}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

        with pytest.raises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

    def test_delete_slice(self):
        idx = timedelta_range(start='1 days', periods=10, freq='D', name='idx')

        # prserve freq
        expected_0_2 = timedelta_range(start='4 days', periods=7, freq='D',
                                       name='idx')
        expected_7_9 = timedelta_range(start='1 days', periods=7, freq='D',
                                       name='idx')

        # reset freq to None
        expected_3_5 = TimedeltaIndex(['1 d', '2 d', '3 d',
                                       '7 d', '8 d', '9 d', '10d'],
                                      freq=None, name='idx')

        cases = {(0, 1, 2): expected_0_2,
                 (7, 8, 9): expected_7_9,
                 (3, 4, 5): expected_3_5}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

            result = idx.delete(slice(n[0], n[-1] + 1))
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

    def test_get_loc(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])

        for method in [None, 'pad', 'backfill', 'nearest']:
            assert idx.get_loc(idx[1], method) == 1
            assert idx.get_loc(idx[1].to_pytimedelta(), method) == 1
            assert idx.get_loc(str(idx[1]), method) == 1

        assert idx.get_loc(idx[1], 'pad',
                           tolerance=Timedelta(0)) == 1
        assert idx.get_loc(idx[1], 'pad',
                           tolerance=np.timedelta64(0, 's')) == 1
        assert idx.get_loc(idx[1], 'pad',
                           tolerance=timedelta(0)) == 1

        with tm.assert_raises_regex(ValueError,
                                    'unit abbreviation w/o a number'):
            idx.get_loc(idx[1], method='nearest', tolerance='foo')

        with pytest.raises(
                ValueError,
                match='tolerance size must match'):
            idx.get_loc(idx[1], method='nearest',
                        tolerance=[Timedelta(0).to_timedelta64(),
                                   Timedelta(0).to_timedelta64()])

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            assert idx.get_loc('1 day 1 hour', method) == loc

        # GH 16909
        assert idx.get_loc(idx[1].to_timedelta64()) == 1

        # GH 16896
        assert idx.get_loc('0 days') == 0

    def test_get_loc_nat(self):
        tidx = TimedeltaIndex(['1 days 01:00:00', 'NaT', '2 days 01:00:00'])

        assert tidx.get_loc(pd.NaT) == 1
        assert tidx.get_loc(None) == 1
        assert tidx.get_loc(float('nan')) == 1
        assert tidx.get_loc(np.nan) == 1

    def test_get_indexer(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])
        tm.assert_numpy_array_equal(idx.get_indexer(idx),
                                    np.array([0, 1, 2], dtype=np.intp))

        target = pd.to_timedelta(['-1 hour', '12 hours', '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 1, 1], dtype=np.intp))

        res = idx.get_indexer(target, 'nearest',
                              tolerance=Timedelta('1 hour'))
        tm.assert_numpy_array_equal(res, np.array([0, -1, 1], dtype=np.intp))
