from warnings import catch_warnings
import pytest
import numpy as np
import pandas as pd
from pandas import (Panel, Series, MultiIndex, DataFrame,
                    Timestamp, Index, date_range)
from pandas.util import testing as tm
from pandas.errors import PerformanceWarning, UnsortedIndexError
from pandas.tests.indexing.common import _mklbl


class TestMultiIndexBasic(object):

    def test_iloc_getitem_multiindex2(self):
        # TODO(wesm): fix this
        pytest.skip('this test was being suppressed, '
                    'needs to be fixed')

        arr = np.random.randn(3, 3)
        df = DataFrame(arr, columns=[[2, 2, 4], [6, 8, 10]],
                       index=[[4, 4, 8], [8, 10, 12]])

        rs = df.iloc[2]
        xp = Series(arr[2], index=df.columns)
        tm.assert_series_equal(rs, xp)

        rs = df.iloc[:, 2]
        xp = Series(arr[:, 2], index=df.index)
        tm.assert_series_equal(rs, xp)

        rs = df.iloc[2, 2]
        xp = df.values[2, 2]
        assert rs == xp

        # for multiple items
        # GH 5528
        rs = df.iloc[[0, 1]]
        xp = df.xs(4, drop_level=False)
        tm.assert_frame_equal(rs, xp)

        tup = zip(*[['a', 'a', 'b', 'b'], ['x', 'y', 'x', 'y']])
        index = MultiIndex.from_tuples(tup)
        df = DataFrame(np.random.randn(4, 4), index=index)
        rs = df.iloc[[2, 3]]
        xp = df.xs('b', drop_level=False)
        tm.assert_frame_equal(rs, xp)

    def test_setitem_multiindex(self):
        with catch_warnings(record=True):

            for index_fn in ('ix', 'loc'):

                def assert_equal(a, b):
                    assert a == b

                def check(target, indexers, value, compare_fn, expected=None):
                    fn = getattr(target, index_fn)
                    fn.__setitem__(indexers, value)
                    result = fn.__getitem__(indexers)
                    if expected is None:
                        expected = value
                    compare_fn(result, expected)
                # GH7190
                index = MultiIndex.from_product([np.arange(0, 100),
                                                 np.arange(0, 80)],
                                                names=['time', 'firm'])
                t, n = 0, 2
                df = DataFrame(np.nan, columns=['A', 'w', 'l', 'a', 'x',
                                                'X', 'd', 'profit'],
                               index=index)
                check(target=df, indexers=((t, n), 'X'), value=0,
                      compare_fn=assert_equal)

                df = DataFrame(-999, columns=['A', 'w', 'l', 'a', 'x',
                                              'X', 'd', 'profit'],
                               index=index)
                check(target=df, indexers=((t, n), 'X'), value=1,
                      compare_fn=assert_equal)

                df = DataFrame(columns=['A', 'w', 'l', 'a', 'x',
                                        'X', 'd', 'profit'],
                               index=index)
                check(target=df, indexers=((t, n), 'X'), value=2,
                      compare_fn=assert_equal)

                # gh-7218: assigning with 0-dim arrays
                df = DataFrame(-999, columns=['A', 'w', 'l', 'a', 'x',
                                              'X', 'd', 'profit'],
                               index=index)
                check(target=df,
                      indexers=((t, n), 'X'),
                      value=np.array(3),
                      compare_fn=assert_equal,
                      expected=3, )

                # GH5206
                df = DataFrame(np.arange(25).reshape(5, 5),
                               columns='A,B,C,D,E'.split(','), dtype=float)
                df['F'] = 99
                row_selection = df['A'] % 2 == 0
                col_selection = ['B', 'C']
                with catch_warnings(record=True):
                    df.ix[row_selection, col_selection] = df['F']
                output = DataFrame(99., index=[0, 2, 4], columns=['B', 'C'])
                with catch_warnings(record=True):
                    tm.assert_frame_equal(df.ix[row_selection, col_selection],
                                          output)
                check(target=df,
                      indexers=(row_selection, col_selection),
                      value=df['F'],
                      compare_fn=tm.assert_frame_equal,
                      expected=output, )

                # GH11372
                idx = MultiIndex.from_product([
                    ['A', 'B', 'C'],
                    date_range('2015-01-01', '2015-04-01', freq='MS')])
                cols = MultiIndex.from_product([
                    ['foo', 'bar'],
                    date_range('2016-01-01', '2016-02-01', freq='MS')])

                df = DataFrame(np.random.random((12, 4)),
                               index=idx, columns=cols)

                subidx = MultiIndex.from_tuples(
                    [('A', Timestamp('2015-01-01')),
                     ('A', Timestamp('2015-02-01'))])
                subcols = MultiIndex.from_tuples(
                    [('foo', Timestamp('2016-01-01')),
                     ('foo', Timestamp('2016-02-01'))])

                vals = DataFrame(np.random.random((2, 2)),
                                 index=subidx, columns=subcols)
                check(target=df,
                      indexers=(subidx, subcols),
                      value=vals,
                      compare_fn=tm.assert_frame_equal, )
                # set all columns
                vals = DataFrame(
                    np.random.random((2, 4)), index=subidx, columns=cols)
                check(target=df,
                      indexers=(subidx, slice(None, None, None)),
                      value=vals,
                      compare_fn=tm.assert_frame_equal, )
                # identity
                copy = df.copy()
                check(target=df, indexers=(df.index, df.columns), value=df,
                      compare_fn=tm.assert_frame_equal, expected=copy)

    def test_loc_getitem_series(self):
        # GH14730
        # passing a series as a key with a MultiIndex
        index = MultiIndex.from_product([[1, 2, 3], ['A', 'B', 'C']])
        x = Series(index=index, data=range(9), dtype=np.float64)
        y = Series([1, 3])
        expected = Series(
            data=[0, 1, 2, 6, 7, 8],
            index=MultiIndex.from_product([[1, 3], ['A', 'B', 'C']]),
            dtype=np.float64)
        result = x.loc[y]
        tm.assert_series_equal(result, expected)

        result = x.loc[[1, 3]]
        tm.assert_series_equal(result, expected)

        # GH15424
        y1 = Series([1, 3], index=[1, 2])
        result = x.loc[y1]
        tm.assert_series_equal(result, expected)

        empty = Series(data=[], dtype=np.float64)
        expected = Series([], index=MultiIndex(
            levels=index.levels, labels=[[], []], dtype=np.float64))
        result = x.loc[empty]
        tm.assert_series_equal(result, expected)

    def test_loc_getitem_array(self):
        # GH15434
        # passing an array as a key with a MultiIndex
        index = MultiIndex.from_product([[1, 2, 3], ['A', 'B', 'C']])
        x = Series(index=index, data=range(9), dtype=np.float64)
        y = np.array([1, 3])
        expected = Series(
            data=[0, 1, 2, 6, 7, 8],
            index=MultiIndex.from_product([[1, 3], ['A', 'B', 'C']]),
            dtype=np.float64)
        result = x.loc[y]
        tm.assert_series_equal(result, expected)

        # empty array:
        empty = np.array([])
        expected = Series([], index=MultiIndex(
            levels=index.levels, labels=[[], []], dtype=np.float64))
        result = x.loc[empty]
        tm.assert_series_equal(result, expected)

        # 0-dim array (scalar):
        scalar = np.int64(1)
        expected = Series(
            data=[0, 1, 2],
            index=['A', 'B', 'C'],
            dtype=np.float64)
        result = x.loc[scalar]
        tm.assert_series_equal(result, expected)

    def test_iloc_getitem_multiindex(self):
        mi_labels = DataFrame(np.random.randn(4, 3),
                              columns=[['i', 'i', 'j'], ['A', 'A', 'B']],
                              index=[['i', 'i', 'j', 'k'],
                                     ['X', 'X', 'Y', 'Y']])

        mi_int = DataFrame(np.random.randn(3, 3),
                           columns=[[2, 2, 4], [6, 8, 10]],
                           index=[[4, 4, 8], [8, 10, 12]])

        # the first row
        rs = mi_int.iloc[0]
        with catch_warnings(record=True):
            xp = mi_int.ix[4].ix[8]
        tm.assert_series_equal(rs, xp, check_names=False)
        assert rs.name == (4, 8)
        assert xp.name == 8

        # 2nd (last) columns
        rs = mi_int.iloc[:, 2]
        with catch_warnings(record=True):
            xp = mi_int.ix[:, 2]
        tm.assert_series_equal(rs, xp)

        # corner column
        rs = mi_int.iloc[2, 2]
        with catch_warnings(record=True):
            xp = mi_int.ix[:, 2].ix[2]
        assert rs == xp

        # this is basically regular indexing
        rs = mi_labels.iloc[2, 2]
        with catch_warnings(record=True):
            xp = mi_labels.ix['j'].ix[:, 'j'].ix[0, 0]
        assert rs == xp

    def test_loc_multiindex(self):

        mi_labels = DataFrame(np.random.randn(3, 3),
                              columns=[['i', 'i', 'j'], ['A', 'A', 'B']],
                              index=[['i', 'i', 'j'], ['X', 'X', 'Y']])

        mi_int = DataFrame(np.random.randn(3, 3),
                           columns=[[2, 2, 4], [6, 8, 10]],
                           index=[[4, 4, 8], [8, 10, 12]])

        # the first row
        rs = mi_labels.loc['i']
        with catch_warnings(record=True):
            xp = mi_labels.ix['i']
        tm.assert_frame_equal(rs, xp)

        # 2nd (last) columns
        rs = mi_labels.loc[:, 'j']
        with catch_warnings(record=True):
            xp = mi_labels.ix[:, 'j']
        tm.assert_frame_equal(rs, xp)

        # corner column
        rs = mi_labels.loc['j'].loc[:, 'j']
        with catch_warnings(record=True):
            xp = mi_labels.ix['j'].ix[:, 'j']
        tm.assert_frame_equal(rs, xp)

        # with a tuple
        rs = mi_labels.loc[('i', 'X')]
        with catch_warnings(record=True):
            xp = mi_labels.ix[('i', 'X')]
        tm.assert_frame_equal(rs, xp)

        rs = mi_int.loc[4]
        with catch_warnings(record=True):
            xp = mi_int.ix[4]
        tm.assert_frame_equal(rs, xp)

    def test_getitem_partial_int(self):
        # GH 12416
        # with single item
        l1 = [10, 20]
        l2 = ['a', 'b']
        df = DataFrame(index=range(2),
                       columns=MultiIndex.from_product([l1, l2]))
        expected = DataFrame(index=range(2),
                             columns=l2)
        result = df[20]
        tm.assert_frame_equal(result, expected)

        # with list
        expected = DataFrame(index=range(2),
                             columns=MultiIndex.from_product([l1[1:], l2]))
        result = df[[20]]
        tm.assert_frame_equal(result, expected)

        # missing item:
        with tm.assert_raises_regex(KeyError, '1'):
            df[1]
        with tm.assert_raises_regex(KeyError, r"'\[1\] not in index'"):
            df[[1]]

    def test_loc_multiindex_indexer_none(self):

        # GH6788
        # multi-index indexer is None (meaning take all)
        attributes = ['Attribute' + str(i) for i in range(1)]
        attribute_values = ['Value' + str(i) for i in range(5)]

        index = MultiIndex.from_product([attributes, attribute_values])
        df = 0.1 * np.random.randn(10, 1 * 5) + 0.5
        df = DataFrame(df, columns=index)
        result = df[attributes]
        tm.assert_frame_equal(result, df)

        # GH 7349
        # loc with a multi-index seems to be doing fallback
        df = DataFrame(np.arange(12).reshape(-1, 1),
                       index=MultiIndex.from_product([[1, 2, 3, 4],
                                                      [1, 2, 3]]))

        expected = df.loc[([1, 2], ), :]
        result = df.loc[[1, 2]]
        tm.assert_frame_equal(result, expected)

    def test_loc_multiindex_incomplete(self):

        # GH 7399
        # incomplete indexers
        s = Series(np.arange(15, dtype='int64'),
                   MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.loc[:, 'a':'c']

        result = s.loc[0:4, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        result = s.loc[:4, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        result = s.loc[0:, 'a':'c']
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected)

        # GH 7400
        # multiindexer gettitem with list of indexers skips wrong element
        s = Series(np.arange(15, dtype='int64'),
                   MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.iloc[[6, 7, 8, 12, 13, 14]]
        result = s.loc[2:4:2, 'a':'c']
        tm.assert_series_equal(result, expected)

    def test_multiindex_perf_warn(self):

        df = DataFrame({'jim': [0, 0, 1, 1],
                        'joe': ['x', 'x', 'z', 'y'],
                        'jolie': np.random.rand(4)}).set_index(['jim', 'joe'])

        with tm.assert_produces_warning(PerformanceWarning,
                                        clear=[pd.core.index]):
            df.loc[(1, 'z')]

        df = df.iloc[[2, 1, 3, 0]]
        with tm.assert_produces_warning(PerformanceWarning):
            df.loc[(0, )]

    def test_series_getitem_multiindex(self):

        # GH 6018
        # series regression getitem with a multi-index

        s = Series([1, 2, 3])
        s.index = MultiIndex.from_tuples([(0, 0), (1, 1), (2, 1)])

        result = s[:, 0]
        expected = Series([1], index=[0])
        tm.assert_series_equal(result, expected)

        result = s.loc[:, 1]
        expected = Series([2, 3], index=[1, 2])
        tm.assert_series_equal(result, expected)

        # xs
        result = s.xs(0, level=0)
        expected = Series([1], index=[0])
        tm.assert_series_equal(result, expected)

        result = s.xs(1, level=1)
        expected = Series([2, 3], index=[1, 2])
        tm.assert_series_equal(result, expected)

        # GH6258
        dt = list(date_range('20130903', periods=3))
        idx = MultiIndex.from_product([list('AB'), dt])
        s = Series([1, 3, 4, 1, 3, 4], index=idx)

        result = s.xs('20130903', level=1)
        expected = Series([1, 1], index=list('AB'))
        tm.assert_series_equal(result, expected)

        # GH5684
        idx = MultiIndex.from_tuples([('a', 'one'), ('a', 'two'), ('b', 'one'),
                                      ('b', 'two')])
        s = Series([1, 2, 3, 4], index=idx)
        s.index.set_names(['L1', 'L2'], inplace=True)
        result = s.xs('one', level='L2')
        expected = Series([1, 3], index=['a', 'b'])
        expected.index.set_names(['L1'], inplace=True)
        tm.assert_series_equal(result, expected)

    def test_xs_multiindex(self):

        # GH2903
        columns = MultiIndex.from_tuples(
            [('a', 'foo'), ('a', 'bar'), ('b', 'hello'),
             ('b', 'world')], names=['lvl0', 'lvl1'])
        df = DataFrame(np.random.randn(4, 4), columns=columns)
        df.sort_index(axis=1, inplace=True)
        result = df.xs('a', level='lvl0', axis=1)
        expected = df.iloc[:, 0:2].loc[:, 'a']
        tm.assert_frame_equal(result, expected)

        result = df.xs('foo', level='lvl1', axis=1)
        expected = df.iloc[:, 1:2].copy()
        expected.columns = expected.columns.droplevel('lvl1')
        tm.assert_frame_equal(result, expected)

    def test_multiindex_setitem(self):

        # GH 3738
        # setting with a multi-index right hand side
        arrays = [np.array(['bar', 'bar', 'baz', 'qux', 'qux', 'bar']),
                  np.array(['one', 'two', 'one', 'one', 'two', 'one']),
                  np.arange(0, 6, 1)]

        df_orig = DataFrame(np.random.randn(6, 3), index=arrays,
                            columns=['A', 'B', 'C']).sort_index()

        expected = df_orig.loc[['bar']] * 2
        df = df_orig.copy()
        df.loc[['bar']] *= 2
        tm.assert_frame_equal(df.loc[['bar']], expected)

        # raise because these have differing levels
        def f():
            df.loc['bar'] *= 2

        pytest.raises(TypeError, f)

        # from SO
        # http://stackoverflow.com/questions/24572040/pandas-access-the-level-of-multiindex-for-inplace-operation
        df_orig = DataFrame.from_dict({'price': {
            ('DE', 'Coal', 'Stock'): 2,
            ('DE', 'Gas', 'Stock'): 4,
            ('DE', 'Elec', 'Demand'): 1,
            ('FR', 'Gas', 'Stock'): 5,
            ('FR', 'Solar', 'SupIm'): 0,
            ('FR', 'Wind', 'SupIm'): 0
        }})
        df_orig.index = MultiIndex.from_tuples(df_orig.index,
                                               names=['Sit', 'Com', 'Type'])

        expected = df_orig.copy()
        expected.iloc[[0, 2, 3]] *= 2

        idx = pd.IndexSlice
        df = df_orig.copy()
        df.loc[idx[:, :, 'Stock'], :] *= 2
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[idx[:, :, 'Stock'], 'price'] *= 2
        tm.assert_frame_equal(df, expected)

    def test_getitem_duplicates_multiindex(self):
        # GH 5725 the 'A' happens to be a valid Timestamp so the doesn't raise
        # the appropriate error, only in PY3 of course!

        index = MultiIndex(levels=[['D', 'B', 'C'],
                                   [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                                   [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        arr = np.random.randn(len(index), 1)
        df = DataFrame(arr, index=index, columns=['val'])
        result = df.val['D']
        expected = Series(arr.ravel()[0:3], name='val', index=Index(
            [26, 37, 57], name='day'))
        tm.assert_series_equal(result, expected)

        def f():
            df.val['A']

        pytest.raises(KeyError, f)

        def f():
            df.val['X']

        pytest.raises(KeyError, f)

        # A is treated as a special Timestamp
        index = MultiIndex(levels=[['A', 'B', 'C'],
                                   [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                                   [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        df = DataFrame(arr, index=index, columns=['val'])
        result = df.val['A']
        expected = Series(arr.ravel()[0:3], name='val', index=Index(
            [26, 37, 57], name='day'))
        tm.assert_series_equal(result, expected)

        def f():
            df.val['X']

        pytest.raises(KeyError, f)

        # GH 7866
        # multi-index slicing with missing indexers
        idx = MultiIndex.from_product([['A', 'B', 'C'],
                                       ['foo', 'bar', 'baz']],
                                      names=['one', 'two'])
        s = Series(np.arange(9, dtype='int64'), index=idx).sort_index()

        exp_idx = MultiIndex.from_product([['A'], ['foo', 'bar', 'baz']],
                                          names=['one', 'two'])
        expected = Series(np.arange(3, dtype='int64'),
                          index=exp_idx).sort_index()

        result = s.loc[['A']]
        tm.assert_series_equal(result, expected)
        result = s.loc[['A', 'D']]
        tm.assert_series_equal(result, expected)

        # not any values found
        pytest.raises(KeyError, lambda: s.loc[['D']])

        # empty ok
        result = s.loc[[]]
        expected = s.iloc[[]]
        tm.assert_series_equal(result, expected)

        idx = pd.IndexSlice
        expected = Series([0, 3, 6], index=MultiIndex.from_product(
            [['A', 'B', 'C'], ['foo']], names=['one', 'two'])).sort_index()

        result = s.loc[idx[:, ['foo']]]
        tm.assert_series_equal(result, expected)
        result = s.loc[idx[:, ['foo', 'bah']]]
        tm.assert_series_equal(result, expected)

        # GH 8737
        # empty indexer
        multi_index = MultiIndex.from_product((['foo', 'bar', 'baz'],
                                               ['alpha', 'beta']))
        df = DataFrame(
            np.random.randn(5, 6), index=range(5), columns=multi_index)
        df = df.sort_index(level=0, axis=1)

        expected = DataFrame(index=range(5),
                             columns=multi_index.reindex([])[0])
        result1 = df.loc[:, ([], slice(None))]
        result2 = df.loc[:, (['foo'], [])]
        tm.assert_frame_equal(result1, expected)
        tm.assert_frame_equal(result2, expected)

        # regression from < 0.14.0
        # GH 7914
        df = DataFrame([[np.mean, np.median], ['mean', 'median']],
                       columns=MultiIndex.from_tuples([('functs', 'mean'),
                                                       ('functs', 'median')]),
                       index=['function', 'name'])
        result = df.loc['function', ('functs', 'mean')]
        assert result == np.mean

    def test_multiindex_assignment(self):

        # GH3777 part 2

        # mixed dtype
        df = DataFrame(np.random.randint(5, 10, size=9).reshape(3, 3),
                       columns=list('abc'),
                       index=[[4, 4, 8], [8, 10, 12]])
        df['d'] = np.nan
        arr = np.array([0., 1.])

        with catch_warnings(record=True):
            df.ix[4, 'd'] = arr
            tm.assert_series_equal(df.ix[4, 'd'],
                                   Series(arr, index=[8, 10], name='d'))

        # single dtype
        df = DataFrame(np.random.randint(5, 10, size=9).reshape(3, 3),
                       columns=list('abc'),
                       index=[[4, 4, 8], [8, 10, 12]])

        with catch_warnings(record=True):
            df.ix[4, 'c'] = arr
            exp = Series(arr, index=[8, 10], name='c', dtype='float64')
            tm.assert_series_equal(df.ix[4, 'c'], exp)

        # scalar ok
        with catch_warnings(record=True):
            df.ix[4, 'c'] = 10
            exp = Series(10, index=[8, 10], name='c', dtype='float64')
            tm.assert_series_equal(df.ix[4, 'c'], exp)

        # invalid assignments
        def f():
            with catch_warnings(record=True):
                df.ix[4, 'c'] = [0, 1, 2, 3]

        pytest.raises(ValueError, f)

        def f():
            with catch_warnings(record=True):
                df.ix[4, 'c'] = [0]

        pytest.raises(ValueError, f)

        # groupby example
        NUM_ROWS = 100
        NUM_COLS = 10
        col_names = ['A' + num for num in
                     map(str, np.arange(NUM_COLS).tolist())]
        index_cols = col_names[:5]

        df = DataFrame(np.random.randint(5, size=(NUM_ROWS, NUM_COLS)),
                       dtype=np.int64, columns=col_names)
        df = df.set_index(index_cols).sort_index()
        grp = df.groupby(level=index_cols[:4])
        df['new_col'] = np.nan

        f_index = np.arange(5)

        def f(name, df2):
            return Series(np.arange(df2.shape[0]),
                          name=df2.index.values[0]).reindex(f_index)

        # TODO(wesm): unused?
        # new_df = pd.concat([f(name, df2) for name, df2 in grp], axis=1).T

        # we are actually operating on a copy here
        # but in this case, that's ok
        for name, df2 in grp:
            new_vals = np.arange(df2.shape[0])
            with catch_warnings(record=True):
                df.ix[name, 'new_col'] = new_vals

    def test_multiindex_label_slicing_with_negative_step(self):
        s = Series(np.arange(20),
                   MultiIndex.from_product([list('abcde'), np.arange(4)]))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            tm.assert_series_equal(s.loc[l_slc], s.iloc[i_slc])
            tm.assert_series_equal(s[l_slc], s.iloc[i_slc])
            with catch_warnings(record=True):
                tm.assert_series_equal(s.ix[l_slc], s.iloc[i_slc])

        assert_slices_equivalent(SLC[::-1], SLC[::-1])

        assert_slices_equivalent(SLC['d'::-1], SLC[15::-1])
        assert_slices_equivalent(SLC[('d', )::-1], SLC[15::-1])

        assert_slices_equivalent(SLC[:'d':-1], SLC[:11:-1])
        assert_slices_equivalent(SLC[:('d', ):-1], SLC[:11:-1])

        assert_slices_equivalent(SLC['d':'b':-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC[('d', ):'b':-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC['d':('b', ):-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC[('d', ):('b', ):-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC['b':'d':-1], SLC[:0])

        assert_slices_equivalent(SLC[('c', 2)::-1], SLC[10::-1])
        assert_slices_equivalent(SLC[:('c', 2):-1], SLC[:9:-1])
        assert_slices_equivalent(SLC[('e', 0):('c', 2):-1], SLC[16:9:-1])

    def test_multiindex_slice_first_level(self):
        # GH 12697
        freq = ['a', 'b', 'c', 'd']
        idx = MultiIndex.from_product([freq, np.arange(500)])
        df = DataFrame(list(range(2000)), index=idx, columns=['Test'])
        df_slice = df.loc[pd.IndexSlice[:, 30:70], :]
        result = df_slice.loc['a']
        expected = DataFrame(list(range(30, 71)),
                             columns=['Test'], index=range(30, 71))
        tm.assert_frame_equal(result, expected)
        result = df_slice.loc['d']
        expected = DataFrame(list(range(1530, 1571)),
                             columns=['Test'], index=range(30, 71))
        tm.assert_frame_equal(result, expected)

    def test_multiindex_symmetric_difference(self):
        # GH 13490
        idx = MultiIndex.from_product([['a', 'b'], ['A', 'B']],
                                      names=['a', 'b'])
        result = idx ^ idx
        assert result.names == idx.names

        idx2 = idx.copy().rename(['A', 'B'])
        result = idx ^ idx2
        assert result.names == [None, None]

    def test_multiindex_contains_dropped(self):
        # GH 19027
        # test that dropped MultiIndex levels are not in the MultiIndex
        # despite continuing to be in the MultiIndex's levels
        idx = MultiIndex.from_product([[1, 2], [3, 4]])
        assert 2 in idx
        idx = idx.drop(2)

        # drop implementation keeps 2 in the levels
        assert 2 in idx.levels[0]
        # but it should no longer be in the index itself
        assert 2 not in idx

        # also applies to strings
        idx = MultiIndex.from_product([['a', 'b'], ['c', 'd']])
        assert 'a' in idx
        idx = idx.drop('a')
        assert 'a' in idx.levels[0]
        assert 'a' not in idx


class TestMultiIndexSlicers(object):

    def test_per_axis_per_level_getitem(self):

        # GH6134
        # example test case
        ix = MultiIndex.from_product([_mklbl('A', 5), _mklbl('B', 7), _mklbl(
            'C', 4), _mklbl('D', 2)])
        df = DataFrame(np.arange(len(ix.get_values())), index=ix)

        result = df.loc[(slice('A1', 'A3'), slice(None), ['C1', 'C3']), :]
        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (a == 'A1' or a == 'A2' or a == 'A3') and (
                               c == 'C1' or c == 'C3')]]
        tm.assert_frame_equal(result, expected)

        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (a == 'A1' or a == 'A2' or a == 'A3') and (
                               c == 'C1' or c == 'C2' or c == 'C3')]]
        result = df.loc[(slice('A1', 'A3'), slice(None), slice('C1', 'C3')), :]
        tm.assert_frame_equal(result, expected)

        # test multi-index slicing with per axis and per index controls
        index = MultiIndex.from_tuples([('A', 1), ('A', 2),
                                        ('A', 3), ('B', 1)],
                                       names=['one', 'two'])
        columns = MultiIndex.from_tuples([('a', 'foo'), ('a', 'bar'),
                                          ('b', 'foo'), ('b', 'bah')],
                                         names=['lvl0', 'lvl1'])

        df = DataFrame(
            np.arange(16, dtype='int64').reshape(
                4, 4), index=index, columns=columns)
        df = df.sort_index(axis=0).sort_index(axis=1)

        # identity
        result = df.loc[(slice(None), slice(None)), :]
        tm.assert_frame_equal(result, df)
        result = df.loc[(slice(None), slice(None)), (slice(None), slice(None))]
        tm.assert_frame_equal(result, df)
        result = df.loc[:, (slice(None), slice(None))]
        tm.assert_frame_equal(result, df)

        # index
        result = df.loc[(slice(None), [1]), :]
        expected = df.iloc[[0, 3]]
        tm.assert_frame_equal(result, expected)

        result = df.loc[(slice(None), 1), :]
        expected = df.iloc[[0, 3]]
        tm.assert_frame_equal(result, expected)

        # columns
        result = df.loc[:, (slice(None), ['foo'])]
        expected = df.iloc[:, [1, 3]]
        tm.assert_frame_equal(result, expected)

        # both
        result = df.loc[(slice(None), 1), (slice(None), ['foo'])]
        expected = df.iloc[[0, 3], [1, 3]]
        tm.assert_frame_equal(result, expected)

        result = df.loc['A', 'a']
        expected = DataFrame(dict(bar=[1, 5, 9], foo=[0, 4, 8]),
                             index=Index([1, 2, 3], name='two'),
                             columns=Index(['bar', 'foo'], name='lvl1'))
        tm.assert_frame_equal(result, expected)

        result = df.loc[(slice(None), [1, 2]), :]
        expected = df.iloc[[0, 1, 3]]
        tm.assert_frame_equal(result, expected)

        # multi-level series
        s = Series(np.arange(len(ix.get_values())), index=ix)
        result = s.loc['A1':'A3', :, ['C1', 'C3']]
        expected = s.loc[[tuple([a, b, c, d])
                          for a, b, c, d in s.index.values
                          if (a == 'A1' or a == 'A2' or a == 'A3') and (
                              c == 'C1' or c == 'C3')]]
        tm.assert_series_equal(result, expected)

        # boolean indexers
        result = df.loc[(slice(None), df.loc[:, ('a', 'bar')] > 5), :]
        expected = df.iloc[[2, 3]]
        tm.assert_frame_equal(result, expected)

        def f():
            df.loc[(slice(None), np.array([True, False])), :]

        pytest.raises(ValueError, f)

        # ambiguous cases
        # these can be multiply interpreted (e.g. in this case
        # as df.loc[slice(None),[1]] as well
        pytest.raises(KeyError, lambda: df.loc[slice(None), [1]])

        result = df.loc[(slice(None), [1]), :]
        expected = df.iloc[[0, 3]]
        tm.assert_frame_equal(result, expected)

        # not lexsorted
        assert df.index.lexsort_depth == 2
        df = df.sort_index(level=1, axis=0)
        assert df.index.lexsort_depth == 0
        with tm.assert_raises_regex(
                UnsortedIndexError,
                'MultiIndex slicing requires the index to be '
                r'lexsorted: slicing on levels \[1\], lexsort depth 0'):
            df.loc[(slice(None), slice('bar')), :]

        # GH 16734: not sorted, but no real slicing
        result = df.loc[(slice(None), df.loc[:, ('a', 'bar')] > 5), :]
        tm.assert_frame_equal(result, df.iloc[[1, 3], :])

    def test_multiindex_slicers_non_unique(self):

        # GH 7106
        # non-unique mi index support
        df = (DataFrame(dict(A=['foo', 'foo', 'foo', 'foo'],
                             B=['a', 'a', 'a', 'a'],
                             C=[1, 2, 1, 3],
                             D=[1, 2, 3, 4]))
              .set_index(['A', 'B', 'C']).sort_index())
        assert not df.index.is_unique
        expected = (DataFrame(dict(A=['foo', 'foo'], B=['a', 'a'],
                                   C=[1, 1], D=[1, 3]))
                    .set_index(['A', 'B', 'C']).sort_index())
        result = df.loc[(slice(None), slice(None), 1), :]
        tm.assert_frame_equal(result, expected)

        # this is equivalent of an xs expression
        result = df.xs(1, level=2, drop_level=False)
        tm.assert_frame_equal(result, expected)

        df = (DataFrame(dict(A=['foo', 'foo', 'foo', 'foo'],
                             B=['a', 'a', 'a', 'a'],
                             C=[1, 2, 1, 2],
                             D=[1, 2, 3, 4]))
              .set_index(['A', 'B', 'C']).sort_index())
        assert not df.index.is_unique
        expected = (DataFrame(dict(A=['foo', 'foo'], B=['a', 'a'],
                                   C=[1, 1], D=[1, 3]))
                    .set_index(['A', 'B', 'C']).sort_index())
        result = df.loc[(slice(None), slice(None), 1), :]
        assert not result.index.is_unique
        tm.assert_frame_equal(result, expected)

        # GH12896
        # numpy-implementation dependent bug
        ints = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 14, 16,
                17, 18, 19, 200000, 200000]
        n = len(ints)
        idx = MultiIndex.from_arrays([['a'] * n, ints])
        result = Series([1] * n, index=idx)
        result = result.sort_index()
        result = result.loc[(slice(None), slice(100000))]
        expected = Series([1] * (n - 2), index=idx[:-2]).sort_index()
        tm.assert_series_equal(result, expected)

    def test_multiindex_slicers_datetimelike(self):

        # GH 7429
        # buggy/inconsistent behavior when slicing with datetime-like
        import datetime
        dates = [datetime.datetime(2012, 1, 1, 12, 12, 12) +
                 datetime.timedelta(days=i) for i in range(6)]
        freq = [1, 2]
        index = MultiIndex.from_product(
            [dates, freq], names=['date', 'frequency'])

        df = DataFrame(
            np.arange(6 * 2 * 4, dtype='int64').reshape(
                -1, 4), index=index, columns=list('ABCD'))

        # multi-axis slicing
        idx = pd.IndexSlice
        expected = df.iloc[[0, 2, 4], [0, 1]]
        result = df.loc[(slice(Timestamp('2012-01-01 12:12:12'),
                               Timestamp('2012-01-03 12:12:12')),
                         slice(1, 1)), slice('A', 'B')]
        tm.assert_frame_equal(result, expected)

        result = df.loc[(idx[Timestamp('2012-01-01 12:12:12'):Timestamp(
            '2012-01-03 12:12:12')], idx[1:1]), slice('A', 'B')]
        tm.assert_frame_equal(result, expected)

        result = df.loc[(slice(Timestamp('2012-01-01 12:12:12'),
                               Timestamp('2012-01-03 12:12:12')), 1),
                        slice('A', 'B')]
        tm.assert_frame_equal(result, expected)

        # with strings
        result = df.loc[(slice('2012-01-01 12:12:12', '2012-01-03 12:12:12'),
                         slice(1, 1)), slice('A', 'B')]
        tm.assert_frame_equal(result, expected)

        result = df.loc[(idx['2012-01-01 12:12:12':'2012-01-03 12:12:12'], 1),
                        idx['A', 'B']]
        tm.assert_frame_equal(result, expected)

    def test_multiindex_slicers_edges(self):
        # GH 8132
        # various edge cases
        df = DataFrame(
            {'A': ['A0'] * 5 + ['A1'] * 5 + ['A2'] * 5,
             'B': ['B0', 'B0', 'B1', 'B1', 'B2'] * 3,
             'DATE': ["2013-06-11", "2013-07-02", "2013-07-09", "2013-07-30",
                      "2013-08-06", "2013-06-11", "2013-07-02", "2013-07-09",
                      "2013-07-30", "2013-08-06", "2013-09-03", "2013-10-01",
                      "2013-07-09", "2013-08-06", "2013-09-03"],
             'VALUES': [22, 35, 14, 9, 4, 40, 18, 4, 2, 5, 1, 2, 3, 4, 2]})

        df['DATE'] = pd.to_datetime(df['DATE'])
        df1 = df.set_index(['A', 'B', 'DATE'])
        df1 = df1.sort_index()

        # A1 - Get all values under "A0" and "A1"
        result = df1.loc[(slice('A1')), :]
        expected = df1.iloc[0:10]
        tm.assert_frame_equal(result, expected)

        # A2 - Get all values from the start to "A2"
        result = df1.loc[(slice('A2')), :]
        expected = df1
        tm.assert_frame_equal(result, expected)

        # A3 - Get all values under "B1" or "B2"
        result = df1.loc[(slice(None), slice('B1', 'B2')), :]
        expected = df1.iloc[[2, 3, 4, 7, 8, 9, 12, 13, 14]]
        tm.assert_frame_equal(result, expected)

        # A4 - Get all values between 2013-07-02 and 2013-07-09
        result = df1.loc[(slice(None), slice(None),
                          slice('20130702', '20130709')), :]
        expected = df1.iloc[[1, 2, 6, 7, 12]]
        tm.assert_frame_equal(result, expected)

        # B1 - Get all values in B0 that are also under A0, A1 and A2
        result = df1.loc[(slice('A2'), slice('B0')), :]
        expected = df1.iloc[[0, 1, 5, 6, 10, 11]]
        tm.assert_frame_equal(result, expected)

        # B2 - Get all values in B0, B1 and B2 (similar to what #2 is doing for
        # the As)
        result = df1.loc[(slice(None), slice('B2')), :]
        expected = df1
        tm.assert_frame_equal(result, expected)

        # B3 - Get all values from B1 to B2 and up to 2013-08-06
        result = df1.loc[(slice(None), slice('B1', 'B2'),
                          slice('2013-08-06')), :]
        expected = df1.iloc[[2, 3, 4, 7, 8, 9, 12, 13]]
        tm.assert_frame_equal(result, expected)

        # B4 - Same as A4 but the start of the date slice is not a key.
        #      shows indexing on a partial selection slice
        result = df1.loc[(slice(None), slice(None),
                          slice('20130701', '20130709')), :]
        expected = df1.iloc[[1, 2, 6, 7, 12]]
        tm.assert_frame_equal(result, expected)

    def test_per_axis_per_level_doc_examples(self):

        # test index maker
        idx = pd.IndexSlice

        # from indexing.rst / advanced
        index = MultiIndex.from_product([_mklbl('A', 4), _mklbl('B', 2),
                                         _mklbl('C', 4), _mklbl('D', 2)])
        columns = MultiIndex.from_tuples([('a', 'foo'), ('a', 'bar'),
                                          ('b', 'foo'), ('b', 'bah')],
                                         names=['lvl0', 'lvl1'])
        df = DataFrame(np.arange(len(index) * len(columns), dtype='int64')
                       .reshape((len(index), len(columns))),
                       index=index, columns=columns)
        result = df.loc[(slice('A1', 'A3'), slice(None), ['C1', 'C3']), :]
        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (a == 'A1' or a == 'A2' or a == 'A3') and (
                               c == 'C1' or c == 'C3')]]
        tm.assert_frame_equal(result, expected)
        result = df.loc[idx['A1':'A3', :, ['C1', 'C3']], :]
        tm.assert_frame_equal(result, expected)

        result = df.loc[(slice(None), slice(None), ['C1', 'C3']), :]
        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (c == 'C1' or c == 'C3')]]
        tm.assert_frame_equal(result, expected)
        result = df.loc[idx[:, :, ['C1', 'C3']], :]
        tm.assert_frame_equal(result, expected)

        # not sorted
        def f():
            df.loc['A1', ('a', slice('foo'))]

        pytest.raises(UnsortedIndexError, f)

        # GH 16734: not sorted, but no real slicing
        tm.assert_frame_equal(df.loc['A1', (slice(None), 'foo')],
                              df.loc['A1'].iloc[:, [0, 2]])

        df = df.sort_index(axis=1)

        # slicing
        df.loc['A1', (slice(None), 'foo')]
        df.loc[(slice(None), slice(None), ['C1', 'C3']), (slice(None), 'foo')]

        # setitem
        df.loc(axis=0)[:, :, ['C1', 'C3']] = -10

    def test_loc_axis_arguments(self):

        index = MultiIndex.from_product([_mklbl('A', 4), _mklbl('B', 2),
                                         _mklbl('C', 4), _mklbl('D', 2)])
        columns = MultiIndex.from_tuples([('a', 'foo'), ('a', 'bar'),
                                          ('b', 'foo'), ('b', 'bah')],
                                         names=['lvl0', 'lvl1'])
        df = DataFrame(np.arange(len(index) * len(columns), dtype='int64')
                       .reshape((len(index), len(columns))),
                       index=index,
                       columns=columns).sort_index().sort_index(axis=1)

        # axis 0
        result = df.loc(axis=0)['A1':'A3', :, ['C1', 'C3']]
        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (a == 'A1' or a == 'A2' or a == 'A3') and (
                               c == 'C1' or c == 'C3')]]
        tm.assert_frame_equal(result, expected)

        result = df.loc(axis='index')[:, :, ['C1', 'C3']]
        expected = df.loc[[tuple([a, b, c, d])
                           for a, b, c, d in df.index.values
                           if (c == 'C1' or c == 'C3')]]
        tm.assert_frame_equal(result, expected)

        # axis 1
        result = df.loc(axis=1)[:, 'foo']
        expected = df.loc[:, (slice(None), 'foo')]
        tm.assert_frame_equal(result, expected)

        result = df.loc(axis='columns')[:, 'foo']
        expected = df.loc[:, (slice(None), 'foo')]
        tm.assert_frame_equal(result, expected)

        # invalid axis
        def f():
            df.loc(axis=-1)[:, :, ['C1', 'C3']]

        pytest.raises(ValueError, f)

        def f():
            df.loc(axis=2)[:, :, ['C1', 'C3']]

        pytest.raises(ValueError, f)

        def f():
            df.loc(axis='foo')[:, :, ['C1', 'C3']]

        pytest.raises(ValueError, f)

    def test_per_axis_per_level_setitem(self):

        # test index maker
        idx = pd.IndexSlice

        # test multi-index slicing with per axis and per index controls
        index = MultiIndex.from_tuples([('A', 1), ('A', 2),
                                        ('A', 3), ('B', 1)],
                                       names=['one', 'two'])
        columns = MultiIndex.from_tuples([('a', 'foo'), ('a', 'bar'),
                                          ('b', 'foo'), ('b', 'bah')],
                                         names=['lvl0', 'lvl1'])

        df_orig = DataFrame(
            np.arange(16, dtype='int64').reshape(
                4, 4), index=index, columns=columns)
        df_orig = df_orig.sort_index(axis=0).sort_index(axis=1)

        # identity
        df = df_orig.copy()
        df.loc[(slice(None), slice(None)), :] = 100
        expected = df_orig.copy()
        expected.iloc[:, :] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc(axis=0)[:, :] = 100
        expected = df_orig.copy()
        expected.iloc[:, :] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None), slice(None)), (slice(None), slice(None))] = 100
        expected = df_orig.copy()
        expected.iloc[:, :] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[:, (slice(None), slice(None))] = 100
        expected = df_orig.copy()
        expected.iloc[:, :] = 100
        tm.assert_frame_equal(df, expected)

        # index
        df = df_orig.copy()
        df.loc[(slice(None), [1]), :] = 100
        expected = df_orig.copy()
        expected.iloc[[0, 3]] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None), 1), :] = 100
        expected = df_orig.copy()
        expected.iloc[[0, 3]] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc(axis=0)[:, 1] = 100
        expected = df_orig.copy()
        expected.iloc[[0, 3]] = 100
        tm.assert_frame_equal(df, expected)

        # columns
        df = df_orig.copy()
        df.loc[:, (slice(None), ['foo'])] = 100
        expected = df_orig.copy()
        expected.iloc[:, [1, 3]] = 100
        tm.assert_frame_equal(df, expected)

        # both
        df = df_orig.copy()
        df.loc[(slice(None), 1), (slice(None), ['foo'])] = 100
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[idx[:, 1], idx[:, ['foo']]] = 100
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] = 100
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc['A', 'a'] = 100
        expected = df_orig.copy()
        expected.iloc[0:3, 0:2] = 100
        tm.assert_frame_equal(df, expected)

        # setting with a list-like
        df = df_orig.copy()
        df.loc[(slice(None), 1), (slice(None), ['foo'])] = np.array(
            [[100, 100], [100, 100]], dtype='int64')
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] = 100
        tm.assert_frame_equal(df, expected)

        # not enough values
        df = df_orig.copy()

        def f():
            df.loc[(slice(None), 1), (slice(None), ['foo'])] = np.array(
                [[100], [100, 100]], dtype='int64')

        pytest.raises(ValueError, f)

        def f():
            df.loc[(slice(None), 1), (slice(None), ['foo'])] = np.array(
                [100, 100, 100, 100], dtype='int64')

        pytest.raises(ValueError, f)

        # with an alignable rhs
        df = df_orig.copy()
        df.loc[(slice(None), 1), (slice(None), ['foo'])] = df.loc[(slice(
            None), 1), (slice(None), ['foo'])] * 5
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] = expected.iloc[[0, 3], [1, 3]] * 5
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None), 1), (slice(None), ['foo'])] *= df.loc[(slice(
            None), 1), (slice(None), ['foo'])]
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] *= expected.iloc[[0, 3], [1, 3]]
        tm.assert_frame_equal(df, expected)

        rhs = df_orig.loc[(slice(None), 1), (slice(None), ['foo'])].copy()
        rhs.loc[:, ('c', 'bah')] = 10
        df = df_orig.copy()
        df.loc[(slice(None), 1), (slice(None), ['foo'])] *= rhs
        expected = df_orig.copy()
        expected.iloc[[0, 3], [1, 3]] *= expected.iloc[[0, 3], [1, 3]]
        tm.assert_frame_equal(df, expected)


class TestMultiIndexPanel(object):

    def test_iloc_getitem_panel_multiindex(self):

        with catch_warnings(record=True):

            # GH 7199
            # Panel with multi-index
            multi_index = MultiIndex.from_tuples([('ONE', 'one'),
                                                  ('TWO', 'two'),
                                                  ('THREE', 'three')],
                                                 names=['UPPER', 'lower'])

            simple_index = [x[0] for x in multi_index]
            wd1 = Panel(items=['First', 'Second'],
                        major_axis=['a', 'b', 'c', 'd'],
                        minor_axis=multi_index)

            wd2 = Panel(items=['First', 'Second'],
                        major_axis=['a', 'b', 'c', 'd'],
                        minor_axis=simple_index)

            expected1 = wd1['First'].iloc[[True, True, True, False], [0, 2]]
            result1 = wd1.iloc[0, [True, True, True, False], [0, 2]]  # WRONG
            tm.assert_frame_equal(result1, expected1)

            expected2 = wd2['First'].iloc[[True, True, True, False], [0, 2]]
            result2 = wd2.iloc[0, [True, True, True, False], [0, 2]]
            tm.assert_frame_equal(result2, expected2)

            expected1 = DataFrame(index=['a'], columns=multi_index,
                                  dtype='float64')
            result1 = wd1.iloc[0, [0], [0, 1, 2]]
            tm.assert_frame_equal(result1, expected1)

            expected2 = DataFrame(index=['a'], columns=simple_index,
                                  dtype='float64')
            result2 = wd2.iloc[0, [0], [0, 1, 2]]
            tm.assert_frame_equal(result2, expected2)

            # GH 7516
            mi = MultiIndex.from_tuples([(0, 'x'), (1, 'y'), (2, 'z')])
            p = Panel(np.arange(3 * 3 * 3, dtype='int64').reshape(3, 3, 3),
                      items=['a', 'b', 'c'], major_axis=mi,
                      minor_axis=['u', 'v', 'w'])
            result = p.iloc[:, 1, 0]
            expected = Series([3, 12, 21], index=['a', 'b', 'c'], name='u')
            tm.assert_series_equal(result, expected)

            result = p.loc[:, (1, 'y'), 'u']
            tm.assert_series_equal(result, expected)

    def test_panel_setitem_with_multiindex(self):

        with catch_warnings(record=True):
            # 10360
            # failing with a multi-index
            arr = np.array([[[1, 2, 3], [0, 0, 0]],
                            [[0, 0, 0], [0, 0, 0]]],
                           dtype=np.float64)

            # reg index
            axes = dict(items=['A', 'B'], major_axis=[0, 1],
                        minor_axis=['X', 'Y', 'Z'])
            p1 = Panel(0., **axes)
            p1.iloc[0, 0, :] = [1, 2, 3]
            expected = Panel(arr, **axes)
            tm.assert_panel_equal(p1, expected)

            # multi-indexes
            axes['items'] = MultiIndex.from_tuples(
                [('A', 'a'), ('B', 'b')])
            p2 = Panel(0., **axes)
            p2.iloc[0, 0, :] = [1, 2, 3]
            expected = Panel(arr, **axes)
            tm.assert_panel_equal(p2, expected)

            axes['major_axis'] = MultiIndex.from_tuples(
                [('A', 1), ('A', 2)])
            p3 = Panel(0., **axes)
            p3.iloc[0, 0, :] = [1, 2, 3]
            expected = Panel(arr, **axes)
            tm.assert_panel_equal(p3, expected)

            axes['minor_axis'] = MultiIndex.from_product(
                [['X'], range(3)])
            p4 = Panel(0., **axes)
            p4.iloc[0, 0, :] = [1, 2, 3]
            expected = Panel(arr, **axes)
            tm.assert_panel_equal(p4, expected)

            arr = np.array(
                [[[1, 0, 0], [2, 0, 0]], [[0, 0, 0], [0, 0, 0]]],
                dtype=np.float64)
            p5 = Panel(0., **axes)
            p5.iloc[0, :, 0] = [1, 2]
            expected = Panel(arr, **axes)
            tm.assert_panel_equal(p5, expected)
