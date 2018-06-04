# pylint: disable=E1103

from warnings import catch_warnings
from numpy.random import randn
import numpy as np
import pytest

import pandas as pd
from pandas.compat import lrange
import pandas.compat as compat
from pandas.util.testing import assert_frame_equal
from pandas import DataFrame, MultiIndex, Series, Index, merge, concat

from pandas._libs import join as libjoin
import pandas.util.testing as tm
from pandas.tests.reshape.merge.test_merge import get_test_data, N, NGROUPS


a_ = np.array


class TestJoin(object):

    def setup_method(self, method):
        # aggregate multiple columns
        self.df = DataFrame({'key1': get_test_data(),
                             'key2': get_test_data(),
                             'data1': np.random.randn(N),
                             'data2': np.random.randn(N)})

        # exclude a couple keys for fun
        self.df = self.df[self.df['key2'] > 1]

        self.df2 = DataFrame({'key1': get_test_data(n=N // 5),
                              'key2': get_test_data(ngroups=NGROUPS // 2,
                                                    n=N // 5),
                              'value': np.random.randn(N // 5)})

        index, data = tm.getMixedTypeDict()
        self.target = DataFrame(data, index=index)

        # Join on string value
        self.source = DataFrame({'MergedA': data['A'], 'MergedD': data['D']},
                                index=data['C'])

    def test_cython_left_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype=np.int64)
        max_group = 5

        ls, rs = libjoin.left_outer_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8, 9, 10])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5, -1, -1])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        tm.assert_numpy_array_equal(ls, exp_ls, check_dtype=False)
        tm.assert_numpy_array_equal(rs, exp_rs, check_dtype=False)

    def test_cython_right_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype=np.int64)
        max_group = 5

        rs, ls = libjoin.left_outer_join(right, left, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        #            0        1        1        1
        exp_li = a_([0, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5,
                     #            2        2        4
                     6, 7, 8, 6, 7, 8, -1])
        exp_ri = a_([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
                     4, 4, 4, 5, 5, 5, 6])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        tm.assert_numpy_array_equal(ls, exp_ls, check_dtype=False)
        tm.assert_numpy_array_equal(rs, exp_rs, check_dtype=False)

    def test_cython_inner_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype=np.int64)
        right = a_([1, 1, 0, 4, 2, 2, 1, 4], dtype=np.int64)
        max_group = 5

        ls, rs = libjoin.inner_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        tm.assert_numpy_array_equal(ls, exp_ls, check_dtype=False)
        tm.assert_numpy_array_equal(rs, exp_rs, check_dtype=False)

    def test_left_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='left')

        joined_both = merge(self.df, self.df2)
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='left')

    def test_right_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='right')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='right')

        joined_both = merge(self.df, self.df2, how='right')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='right')

    def test_full_outer_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='outer')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='outer')

        joined_both = merge(self.df, self.df2, how='outer')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='outer')

    def test_inner_join(self):
        joined_key2 = merge(self.df, self.df2, on='key2', how='inner')
        _check_join(self.df, self.df2, joined_key2, ['key2'], how='inner')

        joined_both = merge(self.df, self.df2, how='inner')
        _check_join(self.df, self.df2, joined_both, ['key1', 'key2'],
                    how='inner')

    def test_handle_overlap(self):
        joined = merge(self.df, self.df2, on='key2',
                       suffixes=['.foo', '.bar'])

        assert 'key1.foo' in joined
        assert 'key1.bar' in joined

    def test_handle_overlap_arbitrary_key(self):
        joined = merge(self.df, self.df2,
                       left_on='key2', right_on='key1',
                       suffixes=['.foo', '.bar'])
        assert 'key1.foo' in joined
        assert 'key2.bar' in joined

    def test_join_on(self):
        target = self.target
        source = self.source

        merged = target.join(source, on='C')
        tm.assert_series_equal(merged['MergedA'], target['A'],
                               check_names=False)
        tm.assert_series_equal(merged['MergedD'], target['D'],
                               check_names=False)

        # join with duplicates (fix regression from DataFrame/Matrix merge)
        df = DataFrame({'key': ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1, 2]}, index=['a', 'b', 'c'])
        joined = df.join(df2, on='key')
        expected = DataFrame({'key': ['a', 'a', 'b', 'b', 'c'],
                              'value': [0, 0, 1, 1, 2]})
        assert_frame_equal(joined, expected)

        # Test when some are missing
        df_a = DataFrame([[1], [2], [3]], index=['a', 'b', 'c'],
                         columns=['one'])
        df_b = DataFrame([['foo'], ['bar']], index=[1, 2],
                         columns=['two'])
        df_c = DataFrame([[1], [2]], index=[1, 2],
                         columns=['three'])
        joined = df_a.join(df_b, on='one')
        joined = joined.join(df_c, on='one')
        assert np.isnan(joined['two']['c'])
        assert np.isnan(joined['three']['c'])

        # merge column not p resent
        pytest.raises(KeyError, target.join, source, on='E')

        # overlap
        source_copy = source.copy()
        source_copy['A'] = 0
        pytest.raises(ValueError, target.join, source_copy, on='A')

    def test_join_on_fails_with_different_right_index(self):
        with pytest.raises(ValueError):
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, left_on='a', right_index=True)

    def test_join_on_fails_with_different_left_index(self):
        with pytest.raises(ValueError):
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)},
                           index=tm.makeCustomIndex(10, 2))
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)})
            merge(df, df2, right_on='b', left_index=True)

    def test_join_on_fails_with_different_column_counts(self):
        with pytest.raises(ValueError):
            df = DataFrame({'a': np.random.choice(['m', 'f'], size=3),
                            'b': np.random.randn(3)})
            df2 = DataFrame({'a': np.random.choice(['m', 'f'], size=10),
                             'b': np.random.randn(10)},
                            index=tm.makeCustomIndex(10, 2))
            merge(df, df2, right_on='a', left_on=['a', 'b'])

    def test_join_on_fails_with_wrong_object_type(self):
        # GH12081
        wrongly_typed = [Series([0, 1]), 2, 'str', None, np.array([0, 1])]
        df = DataFrame({'a': [1, 1]})

        for obj in wrongly_typed:
            with tm.assert_raises_regex(ValueError, str(type(obj))):
                merge(obj, df, left_on='a', right_on='a')
            with tm.assert_raises_regex(ValueError, str(type(obj))):
                merge(df, obj, left_on='a', right_on='a')

    def test_join_on_pass_vector(self):
        expected = self.target.join(self.source, on='C')
        del expected['C']

        join_col = self.target.pop('C')
        result = self.target.join(self.source, on=join_col)
        assert_frame_equal(result, expected)

    def test_join_with_len0(self):
        # nothing to merge
        merged = self.target.join(self.source.reindex([]), on='C')
        for col in self.source:
            assert col in merged
            assert merged[col].isna().all()

        merged2 = self.target.join(self.source.reindex([]), on='C',
                                   how='inner')
        tm.assert_index_equal(merged2.columns, merged.columns)
        assert len(merged2) == 0

    def test_join_on_inner(self):
        df = DataFrame({'key': ['a', 'a', 'd', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1]}, index=['a', 'b'])

        joined = df.join(df2, on='key', how='inner')

        expected = df.join(df2, on='key')
        expected = expected[expected['value'].notna()]
        tm.assert_series_equal(joined['key'], expected['key'],
                               check_dtype=False)
        tm.assert_series_equal(joined['value'], expected['value'],
                               check_dtype=False)
        tm.assert_index_equal(joined.index, expected.index)

    def test_join_on_singlekey_list(self):
        df = DataFrame({'key': ['a', 'a', 'b', 'b', 'c']})
        df2 = DataFrame({'value': [0, 1, 2]}, index=['a', 'b', 'c'])

        # corner cases
        joined = df.join(df2, on=['key'])
        expected = df.join(df2, on='key')

        assert_frame_equal(joined, expected)

    def test_join_on_series(self):
        result = self.target.join(self.source['MergedA'], on='C')
        expected = self.target.join(self.source[['MergedA']], on='C')
        assert_frame_equal(result, expected)

    def test_join_on_series_buglet(self):
        # GH #638
        df = DataFrame({'a': [1, 1]})
        ds = Series([2], index=[1], name='b')
        result = df.join(ds, on='a')
        expected = DataFrame({'a': [1, 1],
                              'b': [2, 2]}, index=df.index)
        tm.assert_frame_equal(result, expected)

    def test_join_index_mixed(self, join_type):
        # no overlapping blocks
        df1 = DataFrame(index=np.arange(10))
        df1['bool'] = True
        df1['string'] = 'foo'

        df2 = DataFrame(index=np.arange(5, 15))
        df2['int'] = 1
        df2['float'] = 1.

        joined = df1.join(df2, how=join_type)
        expected = _join_by_hand(df1, df2, how=join_type)
        assert_frame_equal(joined, expected)

        joined = df2.join(df1, how=join_type)
        expected = _join_by_hand(df2, df1, how=join_type)
        assert_frame_equal(joined, expected)

    def test_join_index_mixed_overlap(self):
        df1 = DataFrame({'A': 1., 'B': 2, 'C': 'foo', 'D': True},
                        index=np.arange(10),
                        columns=['A', 'B', 'C', 'D'])
        assert df1['B'].dtype == np.int64
        assert df1['D'].dtype == np.bool_

        df2 = DataFrame({'A': 1., 'B': 2, 'C': 'foo', 'D': True},
                        index=np.arange(0, 10, 2),
                        columns=['A', 'B', 'C', 'D'])

        # overlap
        joined = df1.join(df2, lsuffix='_one', rsuffix='_two')
        expected_columns = ['A_one', 'B_one', 'C_one', 'D_one',
                            'A_two', 'B_two', 'C_two', 'D_two']
        df1.columns = expected_columns[:4]
        df2.columns = expected_columns[4:]
        expected = _join_by_hand(df1, df2)
        assert_frame_equal(joined, expected)

    def test_join_empty_bug(self):
        # generated an exception in 0.4.3
        x = DataFrame()
        x.join(DataFrame([3], index=[0], columns=['A']), how='outer')

    def test_join_unconsolidated(self):
        # GH #331
        a = DataFrame(randn(30, 2), columns=['a', 'b'])
        c = Series(randn(30))
        a['c'] = c
        d = DataFrame(randn(30, 1), columns=['q'])

        # it works!
        a.join(d)
        d.join(a)

    def test_join_multiindex(self):
        index1 = MultiIndex.from_arrays([['a', 'a', 'a', 'b', 'b', 'b'],
                                         [1, 2, 3, 1, 2, 3]],
                                        names=['first', 'second'])

        index2 = MultiIndex.from_arrays([['b', 'b', 'b', 'c', 'c', 'c'],
                                         [1, 2, 3, 1, 2, 3]],
                                        names=['first', 'second'])

        df1 = DataFrame(data=np.random.randn(6), index=index1,
                        columns=['var X'])
        df2 = DataFrame(data=np.random.randn(6), index=index2,
                        columns=['var Y'])

        df1 = df1.sort_index(level=0)
        df2 = df2.sort_index(level=0)

        joined = df1.join(df2, how='outer')
        ex_index = Index(index1.values).union(Index(index2.values))
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))
        expected.index.names = index1.names
        assert_frame_equal(joined, expected)
        assert joined.index.names == index1.names

        df1 = df1.sort_index(level=1)
        df2 = df2.sort_index(level=1)

        joined = df1.join(df2, how='outer').sort_index(level=0)
        ex_index = Index(index1.values).union(Index(index2.values))
        expected = df1.reindex(ex_index).join(df2.reindex(ex_index))
        expected.index.names = index1.names

        assert_frame_equal(joined, expected)
        assert joined.index.names == index1.names

    def test_join_inner_multiindex(self):
        key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
                'qux', 'snap']
        key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
                'three', 'one']

        data = np.random.randn(len(key1))
        data = DataFrame({'key1': key1, 'key2': key2,
                          'data': data})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        to_join = DataFrame(np.random.randn(10, 3), index=index,
                            columns=['j_one', 'j_two', 'j_three'])

        joined = data.join(to_join, on=['key1', 'key2'], how='inner')
        expected = merge(data, to_join.reset_index(),
                         left_on=['key1', 'key2'],
                         right_on=['first', 'second'], how='inner',
                         sort=False)

        expected2 = merge(to_join, data,
                          right_on=['key1', 'key2'], left_index=True,
                          how='inner', sort=False)
        assert_frame_equal(joined, expected2.reindex_like(joined))

        expected2 = merge(to_join, data, right_on=['key1', 'key2'],
                          left_index=True, how='inner', sort=False)

        expected = expected.drop(['first', 'second'], axis=1)
        expected.index = joined.index

        assert joined.index.is_monotonic
        assert_frame_equal(joined, expected)

        # _assert_same_contents(expected, expected2.loc[:, expected.columns])

    def test_join_hierarchical_mixed(self):
        # GH 2024
        df = DataFrame([(1, 2, 3), (4, 5, 6)], columns=['a', 'b', 'c'])
        new_df = df.groupby(['a']).agg({'b': [np.mean, np.sum]})
        other_df = DataFrame(
            [(1, 2, 3), (7, 10, 6)], columns=['a', 'b', 'd'])
        other_df.set_index('a', inplace=True)
        # GH 9455, 12219
        with tm.assert_produces_warning(UserWarning):
            result = merge(new_df, other_df, left_index=True, right_index=True)
        assert ('b', 'mean') in result
        assert 'b' in result

    def test_join_float64_float32(self):

        a = DataFrame(randn(10, 2), columns=['a', 'b'], dtype=np.float64)
        b = DataFrame(randn(10, 1), columns=['c'], dtype=np.float32)
        joined = a.join(b)
        assert joined.dtypes['a'] == 'float64'
        assert joined.dtypes['b'] == 'float64'
        assert joined.dtypes['c'] == 'float32'

        a = np.random.randint(0, 5, 100).astype('int64')
        b = np.random.random(100).astype('float64')
        c = np.random.random(100).astype('float32')
        df = DataFrame({'a': a, 'b': b, 'c': c})
        xpdf = DataFrame({'a': a, 'b': b, 'c': c})
        s = DataFrame(np.random.random(5).astype('float32'), columns=['md'])
        rs = df.merge(s, left_on='a', right_index=True)
        assert rs.dtypes['a'] == 'int64'
        assert rs.dtypes['b'] == 'float64'
        assert rs.dtypes['c'] == 'float32'
        assert rs.dtypes['md'] == 'float32'

        xp = xpdf.merge(s, left_on='a', right_index=True)
        assert_frame_equal(rs, xp)

    def test_join_many_non_unique_index(self):
        df1 = DataFrame({"a": [1, 1], "b": [1, 1], "c": [10, 20]})
        df2 = DataFrame({"a": [1, 1], "b": [1, 2], "d": [100, 200]})
        df3 = DataFrame({"a": [1, 1], "b": [1, 2], "e": [1000, 2000]})
        idf1 = df1.set_index(["a", "b"])
        idf2 = df2.set_index(["a", "b"])
        idf3 = df3.set_index(["a", "b"])

        result = idf1.join([idf2, idf3], how='outer')

        df_partially_merged = merge(df1, df2, on=['a', 'b'], how='outer')
        expected = merge(df_partially_merged, df3, on=['a', 'b'], how='outer')

        result = result.reset_index()
        expected = expected[result.columns]
        expected['a'] = expected.a.astype('int64')
        expected['b'] = expected.b.astype('int64')
        assert_frame_equal(result, expected)

        df1 = DataFrame({"a": [1, 1, 1], "b": [1, 1, 1], "c": [10, 20, 30]})
        df2 = DataFrame({"a": [1, 1, 1], "b": [1, 1, 2], "d": [100, 200, 300]})
        df3 = DataFrame(
            {"a": [1, 1, 1], "b": [1, 1, 2], "e": [1000, 2000, 3000]})
        idf1 = df1.set_index(["a", "b"])
        idf2 = df2.set_index(["a", "b"])
        idf3 = df3.set_index(["a", "b"])
        result = idf1.join([idf2, idf3], how='inner')

        df_partially_merged = merge(df1, df2, on=['a', 'b'], how='inner')
        expected = merge(df_partially_merged, df3, on=['a', 'b'], how='inner')

        result = result.reset_index()

        assert_frame_equal(result, expected.loc[:, result.columns])

        # GH 11519
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
                        'C': np.random.randn(8),
                        'D': np.random.randn(8)})
        s = Series(np.repeat(np.arange(8), 2),
                   index=np.repeat(np.arange(8), 2), name='TEST')
        inner = df.join(s, how='inner')
        outer = df.join(s, how='outer')
        left = df.join(s, how='left')
        right = df.join(s, how='right')
        assert_frame_equal(inner, outer)
        assert_frame_equal(inner, left)
        assert_frame_equal(inner, right)

    def test_join_sort(self):
        left = DataFrame({'key': ['foo', 'bar', 'baz', 'foo'],
                          'value': [1, 2, 3, 4]})
        right = DataFrame({'value2': ['a', 'b', 'c']},
                          index=['bar', 'baz', 'foo'])

        joined = left.join(right, on='key', sort=True)
        expected = DataFrame({'key': ['bar', 'baz', 'foo', 'foo'],
                              'value': [2, 3, 1, 4],
                              'value2': ['a', 'b', 'c', 'c']},
                             index=[1, 2, 0, 3])
        assert_frame_equal(joined, expected)

        # smoke test
        joined = left.join(right, on='key', sort=False)
        tm.assert_index_equal(joined.index, pd.Index(lrange(4)))

    def test_join_mixed_non_unique_index(self):
        # GH 12814, unorderable types in py3 with a non-unique index
        df1 = DataFrame({'a': [1, 2, 3, 4]}, index=[1, 2, 3, 'a'])
        df2 = DataFrame({'b': [5, 6, 7, 8]}, index=[1, 3, 3, 4])
        result = df1.join(df2)
        expected = DataFrame({'a': [1, 2, 3, 3, 4],
                              'b': [5, np.nan, 6, 7, np.nan]},
                             index=[1, 2, 3, 3, 'a'])
        tm.assert_frame_equal(result, expected)

        df3 = DataFrame({'a': [1, 2, 3, 4]}, index=[1, 2, 2, 'a'])
        df4 = DataFrame({'b': [5, 6, 7, 8]}, index=[1, 2, 3, 4])
        result = df3.join(df4)
        expected = DataFrame({'a': [1, 2, 3, 4], 'b': [5, 6, 6, np.nan]},
                             index=[1, 2, 2, 'a'])
        tm.assert_frame_equal(result, expected)

    def test_join_non_unique_period_index(self):
        # GH #16871
        index = pd.period_range('2016-01-01', periods=16, freq='M')
        df = DataFrame([i for i in range(len(index))],
                       index=index, columns=['pnum'])
        df2 = concat([df, df])
        result = df.join(df2, how='inner', rsuffix='_df2')
        expected = DataFrame(
            np.tile(np.arange(16, dtype=np.int64).repeat(2).reshape(-1, 1), 2),
            columns=['pnum', 'pnum_df2'], index=df2.sort_index().index)
        tm.assert_frame_equal(result, expected)

    def test_mixed_type_join_with_suffix(self):
        # GH #916
        df = DataFrame(np.random.randn(20, 6),
                       columns=['a', 'b', 'c', 'd', 'e', 'f'])
        df.insert(0, 'id', 0)
        df.insert(5, 'dt', 'foo')

        grouped = df.groupby('id')
        mn = grouped.mean()
        cn = grouped.count()

        # it works!
        mn.join(cn, rsuffix='_right')

    def test_join_many(self):
        df = DataFrame(np.random.randn(10, 6), columns=list('abcdef'))
        df_list = [df[['a', 'b']], df[['c', 'd']], df[['e', 'f']]]

        joined = df_list[0].join(df_list[1:])
        tm.assert_frame_equal(joined, df)

        df_list = [df[['a', 'b']][:-2],
                   df[['c', 'd']][2:], df[['e', 'f']][1:9]]

        def _check_diff_index(df_list, result, exp_index):
            reindexed = [x.reindex(exp_index) for x in df_list]
            expected = reindexed[0].join(reindexed[1:])
            tm.assert_frame_equal(result, expected)

        # different join types
        joined = df_list[0].join(df_list[1:], how='outer')
        _check_diff_index(df_list, joined, df.index)

        joined = df_list[0].join(df_list[1:])
        _check_diff_index(df_list, joined, df_list[0].index)

        joined = df_list[0].join(df_list[1:], how='inner')
        _check_diff_index(df_list, joined, df.index[2:8])

        pytest.raises(ValueError, df_list[0].join, df_list[1:], on='a')

    def test_join_many_mixed(self):
        df = DataFrame(np.random.randn(8, 4), columns=['A', 'B', 'C', 'D'])
        df['key'] = ['foo', 'bar'] * 4
        df1 = df.loc[:, ['A', 'B']]
        df2 = df.loc[:, ['C', 'D']]
        df3 = df.loc[:, ['key']]

        result = df1.join([df2, df3])
        assert_frame_equal(result, df)

    def test_join_dups(self):

        # joining dups
        df = concat([DataFrame(np.random.randn(10, 4),
                               columns=['A', 'A', 'B', 'B']),
                     DataFrame(np.random.randint(0, 10, size=20)
                               .reshape(10, 2),
                               columns=['A', 'C'])],
                    axis=1)

        expected = concat([df, df], axis=1)
        result = df.join(df, rsuffix='_2')
        result.columns = expected.columns
        assert_frame_equal(result, expected)

        # GH 4975, invalid join on dups
        w = DataFrame(np.random.randn(4, 2), columns=["x", "y"])
        x = DataFrame(np.random.randn(4, 2), columns=["x", "y"])
        y = DataFrame(np.random.randn(4, 2), columns=["x", "y"])
        z = DataFrame(np.random.randn(4, 2), columns=["x", "y"])

        dta = x.merge(y, left_index=True, right_index=True).merge(
            z, left_index=True, right_index=True, how="outer")
        dta = dta.merge(w, left_index=True, right_index=True)
        expected = concat([x, y, z, w], axis=1)
        expected.columns = ['x_x', 'y_x', 'x_y',
                            'y_y', 'x_x', 'y_x', 'x_y', 'y_y']
        assert_frame_equal(dta, expected)

    def test_panel_join(self):
        with catch_warnings(record=True):
            panel = tm.makePanel()
            tm.add_nans(panel)

            p1 = panel.iloc[:2, :10, :3]
            p2 = panel.iloc[2:, 5:, 2:]

            # left join
            result = p1.join(p2)
            expected = p1.copy()
            expected['ItemC'] = p2['ItemC']
            tm.assert_panel_equal(result, expected)

            # right join
            result = p1.join(p2, how='right')
            expected = p2.copy()
            expected['ItemA'] = p1['ItemA']
            expected['ItemB'] = p1['ItemB']
            expected = expected.reindex(items=['ItemA', 'ItemB', 'ItemC'])
            tm.assert_panel_equal(result, expected)

            # inner join
            result = p1.join(p2, how='inner')
            expected = panel.iloc[:, 5:10, 2:3]
            tm.assert_panel_equal(result, expected)

            # outer join
            result = p1.join(p2, how='outer')
            expected = p1.reindex(major=panel.major_axis,
                                  minor=panel.minor_axis)
            expected = expected.join(p2.reindex(major=panel.major_axis,
                                                minor=panel.minor_axis))
            tm.assert_panel_equal(result, expected)

    def test_panel_join_overlap(self):
        with catch_warnings(record=True):
            panel = tm.makePanel()
            tm.add_nans(panel)

            p1 = panel.loc[['ItemA', 'ItemB', 'ItemC']]
            p2 = panel.loc[['ItemB', 'ItemC']]

            # Expected index is
            #
            # ItemA, ItemB_p1, ItemC_p1, ItemB_p2, ItemC_p2
            joined = p1.join(p2, lsuffix='_p1', rsuffix='_p2')
            p1_suf = p1.loc[['ItemB', 'ItemC']].add_suffix('_p1')
            p2_suf = p2.loc[['ItemB', 'ItemC']].add_suffix('_p2')
            no_overlap = panel.loc[['ItemA']]
            expected = no_overlap.join(p1_suf.join(p2_suf))
            tm.assert_panel_equal(joined, expected)

    def test_panel_join_many(self):
        with catch_warnings(record=True):
            tm.K = 10
            panel = tm.makePanel()
            tm.K = 4

            panels = [panel.iloc[:2], panel.iloc[2:6], panel.iloc[6:]]

            joined = panels[0].join(panels[1:])
            tm.assert_panel_equal(joined, panel)

            panels = [panel.iloc[:2, :-5],
                      panel.iloc[2:6, 2:],
                      panel.iloc[6:, 5:-7]]

            data_dict = {}
            for p in panels:
                data_dict.update(p.iteritems())

            joined = panels[0].join(panels[1:], how='inner')
            expected = pd.Panel.from_dict(data_dict, intersect=True)
            tm.assert_panel_equal(joined, expected)

            joined = panels[0].join(panels[1:], how='outer')
            expected = pd.Panel.from_dict(data_dict, intersect=False)
            tm.assert_panel_equal(joined, expected)

            # edge cases
            pytest.raises(ValueError, panels[0].join, panels[1:],
                          how='outer', lsuffix='foo', rsuffix='bar')
            pytest.raises(ValueError, panels[0].join, panels[1:],
                          how='right')


def _check_join(left, right, result, join_col, how='left',
                lsuffix='_x', rsuffix='_y'):

    # some smoke tests
    for c in join_col:
        assert(result[c].notna().all())

    left_grouped = left.groupby(join_col)
    right_grouped = right.groupby(join_col)

    for group_key, group in result.groupby(join_col):
        l_joined = _restrict_to_columns(group, left.columns, lsuffix)
        r_joined = _restrict_to_columns(group, right.columns, rsuffix)

        try:
            lgroup = left_grouped.get_group(group_key)
        except KeyError:
            if how in ('left', 'inner'):
                raise AssertionError('key %s should not have been in the join'
                                     % str(group_key))

            _assert_all_na(l_joined, left.columns, join_col)
        else:
            _assert_same_contents(l_joined, lgroup)

        try:
            rgroup = right_grouped.get_group(group_key)
        except KeyError:
            if how in ('right', 'inner'):
                raise AssertionError('key %s should not have been in the join'
                                     % str(group_key))

            _assert_all_na(r_joined, right.columns, join_col)
        else:
            _assert_same_contents(r_joined, rgroup)


def _restrict_to_columns(group, columns, suffix):
    found = [c for c in group.columns
             if c in columns or c.replace(suffix, '') in columns]

    # filter
    group = group.loc[:, found]

    # get rid of suffixes, if any
    group = group.rename(columns=lambda x: x.replace(suffix, ''))

    # put in the right order...
    group = group.loc[:, columns]

    return group


def _assert_same_contents(join_chunk, source):
    NA_SENTINEL = -1234567  # drop_duplicates not so NA-friendly...

    jvalues = join_chunk.fillna(NA_SENTINEL).drop_duplicates().values
    svalues = source.fillna(NA_SENTINEL).drop_duplicates().values

    rows = {tuple(row) for row in jvalues}
    assert(len(rows) == len(source))
    assert(all(tuple(row) in rows for row in svalues))


def _assert_all_na(join_chunk, source_columns, join_col):
    for c in source_columns:
        if c in join_col:
            continue
        assert(join_chunk[c].isna().all())


def _join_by_hand(a, b, how='left'):
    join_index = a.index.join(b.index, how=how)

    a_re = a.reindex(join_index)
    b_re = b.reindex(join_index)

    result_columns = a.columns.append(b.columns)

    for col, s in compat.iteritems(b_re):
        a_re[col] = s
    return a_re.reindex(columns=result_columns)
