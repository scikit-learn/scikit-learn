""" test indexing with ix """

import pytest

from warnings import catch_warnings

import numpy as np
import pandas as pd

from pandas.core.dtypes.common import is_scalar
from pandas.compat import lrange
from pandas import Series, DataFrame, option_context, MultiIndex
from pandas.util import testing as tm
from pandas.errors import PerformanceWarning


class TestIX(object):

    def test_ix_deprecation(self):
        # GH 15114

        df = DataFrame({'A': [1, 2, 3]})
        with tm.assert_produces_warning(DeprecationWarning,
                                        check_stacklevel=False):
            df.ix[1, 'A']

    def test_ix_loc_setitem_consistency(self):

        # GH 5771
        # loc with slice and series
        s = Series(0, index=[4, 5, 6])
        s.loc[4:5] += 1
        expected = Series([1, 1, 0], index=[4, 5, 6])
        tm.assert_series_equal(s, expected)

        # GH 5928
        # chained indexing assignment
        df = DataFrame({'a': [0, 1, 2]})
        expected = df.copy()
        with catch_warnings(record=True):
            expected.ix[[0, 1, 2], 'a'] = -expected.ix[[0, 1, 2], 'a']

        with catch_warnings(record=True):
            df['a'].ix[[0, 1, 2]] = -df['a'].ix[[0, 1, 2]]
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'a': [0, 1, 2], 'b': [0, 1, 2]})
        with catch_warnings(record=True):
            df['a'].ix[[0, 1, 2]] = -df['a'].ix[[0, 1, 2]].astype(
                'float64') + 0.5
        expected = DataFrame({'a': [0.5, -0.5, -1.5], 'b': [0, 1, 2]})
        tm.assert_frame_equal(df, expected)

        # GH 8607
        # ix setitem consistency
        df = DataFrame({'delta': [1174, 904, 161],
                        'elapsed': [7673, 9277, 1470],
                        'timestamp': [1413840976, 1413842580, 1413760580]})
        expected = DataFrame({'delta': [1174, 904, 161],
                              'elapsed': [7673, 9277, 1470],
                              'timestamp': pd.to_datetime(
                                  [1413840976, 1413842580, 1413760580],
                                  unit='s')
                              })

        df2 = df.copy()
        df2['timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

        df2 = df.copy()
        df2.loc[:, 'timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

        df2 = df.copy()
        with catch_warnings(record=True):
            df2.ix[:, 2] = pd.to_datetime(df['timestamp'], unit='s')
        tm.assert_frame_equal(df2, expected)

    def test_ix_loc_consistency(self):

        # GH 8613
        # some edge cases where ix/loc should return the same
        # this is not an exhaustive case

        def compare(result, expected):
            if is_scalar(expected):
                assert result == expected
            else:
                assert expected.equals(result)

        # failure cases for .loc, but these work for .ix
        df = DataFrame(np.random.randn(5, 4), columns=list('ABCD'))
        for key in [slice(1, 3), tuple([slice(0, 2), slice(0, 2)]),
                    tuple([slice(0, 2), df.columns[0:2]])]:

            for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                          tm.makeDateIndex, tm.makePeriodIndex,
                          tm.makeTimedeltaIndex]:
                df.index = index(len(df.index))
                with catch_warnings(record=True):
                    df.ix[key]

                pytest.raises(TypeError, lambda: df.loc[key])

        df = DataFrame(np.random.randn(5, 4), columns=list('ABCD'),
                       index=pd.date_range('2012-01-01', periods=5))

        for key in ['2012-01-03',
                    '2012-01-31',
                    slice('2012-01-03', '2012-01-03'),
                    slice('2012-01-03', '2012-01-04'),
                    slice('2012-01-03', '2012-01-06', 2),
                    slice('2012-01-03', '2012-01-31'),
                    tuple([[True, True, True, False, True]]), ]:

            # getitem

            # if the expected raises, then compare the exceptions
            try:
                with catch_warnings(record=True):
                    expected = df.ix[key]
            except KeyError:
                pytest.raises(KeyError, lambda: df.loc[key])
                continue

            result = df.loc[key]
            compare(result, expected)

            # setitem
            df1 = df.copy()
            df2 = df.copy()

            with catch_warnings(record=True):
                df1.ix[key] = 10
            df2.loc[key] = 10
            compare(df2, df1)

        # edge cases
        s = Series([1, 2, 3, 4], index=list('abde'))

        result1 = s['a':'c']
        with catch_warnings(record=True):
            result2 = s.ix['a':'c']
        result3 = s.loc['a':'c']
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

        # now work rather than raising KeyError
        s = Series(range(5), [-2, -1, 1, 2, 3])

        with catch_warnings(record=True):
            result1 = s.ix[-10:3]
        result2 = s.loc[-10:3]
        tm.assert_series_equal(result1, result2)

        with catch_warnings(record=True):
            result1 = s.ix[0:3]
        result2 = s.loc[0:3]
        tm.assert_series_equal(result1, result2)

    def test_ix_weird_slicing(self):
        # http://stackoverflow.com/q/17056560/1240268
        df = DataFrame({'one': [1, 2, 3, np.nan, np.nan],
                        'two': [1, 2, 3, 4, 5]})
        df.loc[df['one'] > 1, 'two'] = -df['two']

        expected = DataFrame({'one': {0: 1.0,
                                      1: 2.0,
                                      2: 3.0,
                                      3: np.nan,
                                      4: np.nan},
                              'two': {0: 1,
                                      1: -2,
                                      2: -3,
                                      3: 4,
                                      4: 5}})
        tm.assert_frame_equal(df, expected)

    def test_ix_general(self):

        # ix general issues

        # GH 2817
        data = {'amount': {0: 700, 1: 600, 2: 222, 3: 333, 4: 444},
                'col': {0: 3.5, 1: 3.5, 2: 4.0, 3: 4.0, 4: 4.0},
                'year': {0: 2012, 1: 2011, 2: 2012, 3: 2012, 4: 2012}}
        df = DataFrame(data).set_index(keys=['col', 'year'])
        key = 4.0, 2012

        # emits a PerformanceWarning, ok
        with tm.assert_produces_warning(PerformanceWarning):
            tm.assert_frame_equal(df.loc[key], df.iloc[2:])

        # this is ok
        df.sort_index(inplace=True)
        res = df.loc[key]

        # col has float dtype, result should be Float64Index
        index = MultiIndex.from_arrays([[4.] * 3, [2012] * 3],
                                       names=['col', 'year'])
        expected = DataFrame({'amount': [222, 333, 444]}, index=index)
        tm.assert_frame_equal(res, expected)

    def test_ix_assign_column_mixed(self):
        # GH #1142
        df = DataFrame(tm.getSeriesData())
        df['foo'] = 'bar'

        orig = df.loc[:, 'B'].copy()
        df.loc[:, 'B'] = df.loc[:, 'B'] + 1
        tm.assert_series_equal(df.B, orig + 1)

        # GH 3668, mixed frame with series value
        df = DataFrame({'x': lrange(10), 'y': lrange(10, 20), 'z': 'bar'})
        expected = df.copy()

        for i in range(5):
            indexer = i * 2
            v = 1000 + i * 200
            expected.loc[indexer, 'y'] = v
            assert expected.loc[indexer, 'y'] == v

        df.loc[df.x % 2 == 0, 'y'] = df.loc[df.x % 2 == 0, 'y'] * 100
        tm.assert_frame_equal(df, expected)

        # GH 4508, making sure consistency of assignments
        df = DataFrame({'a': [1, 2, 3], 'b': [0, 1, 2]})
        df.loc[[0, 2, ], 'b'] = [100, -100]
        expected = DataFrame({'a': [1, 2, 3], 'b': [100, 1, -100]})
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'a': lrange(4)})
        df['b'] = np.nan
        df.loc[[1, 3], 'b'] = [100, -100]
        expected = DataFrame({'a': [0, 1, 2, 3],
                              'b': [np.nan, 100, np.nan, -100]})
        tm.assert_frame_equal(df, expected)

        # ok, but chained assignments are dangerous
        # if we turn off chained assignment it will work
        with option_context('chained_assignment', None):
            df = DataFrame({'a': lrange(4)})
            df['b'] = np.nan
            df['b'].loc[[1, 3]] = [100, -100]
            tm.assert_frame_equal(df, expected)

    def test_ix_get_set_consistency(self):

        # GH 4544
        # ix/loc get/set not consistent when
        # a mixed int/string index
        df = DataFrame(np.arange(16).reshape((4, 4)),
                       columns=['a', 'b', 8, 'c'],
                       index=['e', 7, 'f', 'g'])

        with catch_warnings(record=True):
            assert df.ix['e', 8] == 2
        assert df.loc['e', 8] == 2

        with catch_warnings(record=True):
            df.ix['e', 8] = 42
            assert df.ix['e', 8] == 42
        assert df.loc['e', 8] == 42

        df.loc['e', 8] = 45
        with catch_warnings(record=True):
            assert df.ix['e', 8] == 45
        assert df.loc['e', 8] == 45

    def test_ix_slicing_strings(self):
        # see gh-3836
        data = {'Classification':
                ['SA EQUITY CFD', 'bbb', 'SA EQUITY', 'SA SSF', 'aaa'],
                'Random': [1, 2, 3, 4, 5],
                'X': ['correct', 'wrong', 'correct', 'correct', 'wrong']}
        df = DataFrame(data)
        x = df[~df.Classification.isin(['SA EQUITY CFD', 'SA EQUITY', 'SA SSF'
                                        ])]
        with catch_warnings(record=True):
            df.ix[x.index, 'X'] = df['Classification']

        expected = DataFrame({'Classification': {0: 'SA EQUITY CFD',
                                                 1: 'bbb',
                                                 2: 'SA EQUITY',
                                                 3: 'SA SSF',
                                                 4: 'aaa'},
                              'Random': {0: 1,
                                         1: 2,
                                         2: 3,
                                         3: 4,
                                         4: 5},
                              'X': {0: 'correct',
                                    1: 'bbb',
                                    2: 'correct',
                                    3: 'correct',
                                    4: 'aaa'}})  # bug was 4: 'bbb'

        tm.assert_frame_equal(df, expected)

    def test_ix_setitem_out_of_bounds_axis_0(self):
        df = DataFrame(
            np.random.randn(2, 5), index=["row%s" % i for i in range(2)],
            columns=["col%s" % i for i in range(5)])
        with catch_warnings(record=True):
            pytest.raises(ValueError, df.ix.__setitem__, (2, 0), 100)

    def test_ix_setitem_out_of_bounds_axis_1(self):
        df = DataFrame(
            np.random.randn(5, 2), index=["row%s" % i for i in range(5)],
            columns=["col%s" % i for i in range(2)])
        with catch_warnings(record=True):
            pytest.raises(ValueError, df.ix.__setitem__, (0, 2), 100)

    def test_ix_empty_list_indexer_is_ok(self):
        with catch_warnings(record=True):
            from pandas.util.testing import makeCustomDataframe as mkdf
            df = mkdf(5, 2)
            # vertical empty
            tm.assert_frame_equal(df.ix[:, []], df.iloc[:, :0],
                                  check_index_type=True,
                                  check_column_type=True)
            # horizontal empty
            tm.assert_frame_equal(df.ix[[], :], df.iloc[:0, :],
                                  check_index_type=True,
                                  check_column_type=True)
            # horizontal empty
            tm.assert_frame_equal(df.ix[[]], df.iloc[:0, :],
                                  check_index_type=True,
                                  check_column_type=True)

    def test_ix_duplicate_returns_series(self):
        df = DataFrame(np.random.randn(3, 3), index=[0.1, 0.2, 0.2],
                       columns=list('abc'))
        with catch_warnings(record=True):
            r = df.ix[0.2, 'a']
        e = df.loc[0.2, 'a']
        tm.assert_series_equal(r, e)
