# -*- coding: utf-8 -*-

import pytest

import pandas as pd
import pandas.compat as compat
import numpy as np
from pandas import (Series, DataFrame, Timestamp, Categorical,
                    CategoricalIndex, Interval, Index)
from pandas.util.testing import assert_series_equal, assert_frame_equal
from pandas.util import testing as tm
from pandas.core.dtypes.common import is_categorical_dtype
from pandas.api.types import CategoricalDtype as CDT
from pandas.core.dtypes.dtypes import CategoricalDtype


class TestCategoricalIndex(object):

    def setup_method(self, method):

        self.df = DataFrame({'A': np.arange(6, dtype='int64'),
                             'B': Series(list('aabbca')).astype(
                                 CDT(list('cab')))}).set_index('B')
        self.df2 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': Series(list('aabbca')).astype(
                                  CDT(list('cabe')))}).set_index('B')
        self.df3 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': (Series([1, 1, 2, 1, 3, 2])
                                    .astype(CDT([3, 2, 1], ordered=True)))
                              }).set_index('B')
        self.df4 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': (Series([1, 1, 2, 1, 3, 2])
                                    .astype(CDT([3, 2, 1], ordered=False)))
                              }).set_index('B')

    def test_loc_scalar(self):
        result = self.df.loc['a']
        expected = (DataFrame({'A': [0, 1, 5],
                               'B': (Series(list('aaa'))
                                     .astype(CDT(list('cab'))))})
                    .set_index('B'))
        assert_frame_equal(result, expected)

        df = self.df.copy()
        df.loc['a'] = 20
        expected = (DataFrame({'A': [20, 20, 2, 3, 4, 20],
                               'B': (Series(list('aabbca'))
                                     .astype(CDT(list('cab'))))})
                    .set_index('B'))
        assert_frame_equal(df, expected)

        # value not in the categories
        pytest.raises(KeyError, lambda: df.loc['d'])

        def f():
            df.loc['d'] = 10

        pytest.raises(TypeError, f)

        def f():
            df.loc['d', 'A'] = 10

        pytest.raises(TypeError, f)

        def f():
            df.loc['d', 'C'] = 10

        pytest.raises(TypeError, f)

    def test_getitem_scalar(self):

        cats = Categorical([Timestamp('12-31-1999'),
                            Timestamp('12-31-2000')])

        s = Series([1, 2], index=cats)

        expected = s.iloc[0]
        result = s[cats[0]]
        assert result == expected

    def test_slicing_directly(self):
        cat = Categorical(["a", "b", "c", "d", "a", "b", "c"])
        sliced = cat[3]
        assert sliced == "d"
        sliced = cat[3:5]
        expected = Categorical(["d", "a"], categories=['a', 'b', 'c', 'd'])
        tm.assert_numpy_array_equal(sliced._codes, expected._codes)
        tm.assert_index_equal(sliced.categories, expected.categories)

    def test_slicing(self):
        cat = Series(Categorical([1, 2, 3, 4]))
        reversed = cat[::-1]
        exp = np.array([4, 3, 2, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(reversed.__array__(), exp)

        df = DataFrame({'value': (np.arange(100) + 1).astype('int64')})
        df['D'] = pd.cut(df.value, bins=[0, 25, 50, 75, 100])

        expected = Series([11, Interval(0, 25)], index=['value', 'D'], name=10)
        result = df.iloc[10]
        tm.assert_series_equal(result, expected)

        expected = DataFrame({'value': np.arange(11, 21).astype('int64')},
                             index=np.arange(10, 20).astype('int64'))
        expected['D'] = pd.cut(expected.value, bins=[0, 25, 50, 75, 100])
        result = df.iloc[10:20]
        tm.assert_frame_equal(result, expected)

        expected = Series([9, Interval(0, 25)], index=['value', 'D'], name=8)
        result = df.loc[8]
        tm.assert_series_equal(result, expected)

    def test_slicing_and_getting_ops(self):

        # systematically test the slicing operations:
        #  for all slicing ops:
        #   - returning a dataframe
        #   - returning a column
        #   - returning a row
        #   - returning a single value

        cats = Categorical(
            ["a", "c", "b", "c", "c", "c", "c"], categories=["a", "b", "c"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n"])
        values = [1, 2, 3, 4, 5, 6, 7]
        df = DataFrame({"cats": cats, "values": values}, index=idx)

        # the expected values
        cats2 = Categorical(["b", "c"], categories=["a", "b", "c"])
        idx2 = Index(["j", "k"])
        values2 = [3, 4]

        # 2:4,: | "j":"k",:
        exp_df = DataFrame({"cats": cats2, "values": values2}, index=idx2)

        # :,"cats" | :,0
        exp_col = Series(cats, index=idx, name='cats')

        # "j",: | 2,:
        exp_row = Series(["b", 3], index=["cats", "values"], dtype="object",
                         name="j")

        # "j","cats | 2,0
        exp_val = "b"

        # iloc
        # frame
        res_df = df.iloc[2:4, :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.iloc[2, :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.iloc[:, 0]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.iloc[2, 0]
        assert res_val == exp_val

        # loc
        # frame
        res_df = df.loc["j":"k", :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.loc["j", :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.loc[:, "cats"]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.loc["j", "cats"]
        assert res_val == exp_val

        # ix
        # frame
        # res_df = df.loc["j":"k",[0,1]] # doesn't work?
        res_df = df.loc["j":"k", :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.loc["j", :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.loc[:, "cats"]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.loc["j", df.columns[0]]
        assert res_val == exp_val

        # iat
        res_val = df.iat[2, 0]
        assert res_val == exp_val

        # at
        res_val = df.at["j", "cats"]
        assert res_val == exp_val

        # fancy indexing
        exp_fancy = df.iloc[[2]]

        res_fancy = df[df["cats"] == "b"]
        tm.assert_frame_equal(res_fancy, exp_fancy)
        res_fancy = df[df["values"] == 3]
        tm.assert_frame_equal(res_fancy, exp_fancy)

        # get_value
        res_val = df.at["j", "cats"]
        assert res_val == exp_val

        # i : int, slice, or sequence of integers
        res_row = df.iloc[2]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        res_df = df.iloc[slice(2, 4)]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        res_df = df.iloc[[2, 3]]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        res_col = df.iloc[:, 0]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        res_df = df.iloc[:, slice(0, 2)]
        tm.assert_frame_equal(res_df, df)
        assert is_categorical_dtype(res_df["cats"])

        res_df = df.iloc[:, [0, 1]]
        tm.assert_frame_equal(res_df, df)
        assert is_categorical_dtype(res_df["cats"])

    def test_slicing_doc_examples(self):

        # GH 7918
        cats = Categorical(["a", "b", "b", "b", "c", "c", "c"],
                           categories=["a", "b", "c"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n", ])
        values = [1, 2, 2, 2, 3, 4, 5]
        df = DataFrame({"cats": cats, "values": values}, index=idx)

        result = df.iloc[2:4, :]
        expected = DataFrame(
            {"cats": Categorical(['b', 'b'], categories=['a', 'b', 'c']),
             "values": [2, 2]}, index=['j', 'k'])
        tm.assert_frame_equal(result, expected)

        result = df.iloc[2:4, :].dtypes
        expected = Series(['category', 'int64'], ['cats', 'values'])
        tm.assert_series_equal(result, expected)

        result = df.loc["h":"j", "cats"]
        expected = Series(Categorical(['a', 'b', 'b'],
                                      categories=['a', 'b', 'c']),
                          index=['h', 'i', 'j'], name='cats')
        tm.assert_series_equal(result, expected)

        result = df.loc["h":"j", df.columns[0:1]]
        expected = DataFrame({'cats': Categorical(['a', 'b', 'b'],
                                                  categories=['a', 'b', 'c'])},
                             index=['h', 'i', 'j'])
        tm.assert_frame_equal(result, expected)

    def test_getitem_category_type(self):
        # GH 14580
        # test iloc() on Series with Categorical data

        s = Series([1, 2, 3]).astype('category')

        # get slice
        result = s.iloc[0:2]
        expected = Series([1, 2]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

        # get list of indexes
        result = s.iloc[[0, 1]]
        expected = Series([1, 2]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

        # get boolean array
        result = s.iloc[[True, False, False]]
        expected = Series([1]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

    def test_loc_listlike(self):

        # list of labels
        result = self.df.loc[['c', 'a']]
        expected = self.df.iloc[[4, 0, 1, 5]]
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.loc[['a', 'b', 'e']]
        exp_index = CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan]}, index=exp_index)
        assert_frame_equal(result, expected, check_index_type=True)

        # element in the categories but not in the values
        pytest.raises(KeyError, lambda: self.df2.loc['e'])

        # assign is ok
        df = self.df2.copy()
        df.loc['e'] = 20
        result = df.loc[['a', 'b', 'e']]
        exp_index = CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, 20]}, index=exp_index)
        assert_frame_equal(result, expected)

        df = self.df2.copy()
        result = df.loc[['a', 'b', 'e']]
        exp_index = CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan]}, index=exp_index)
        assert_frame_equal(result, expected, check_index_type=True)

        # not all labels in the categories
        with pytest.raises(KeyError):
            self.df2.loc[['a', 'd']]

    def test_loc_listlike_dtypes(self):
        # GH 11586

        # unique categories and codes
        index = CategoricalIndex(['a', 'b', 'c'])
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}, index=index)

        # unique slice
        res = df.loc[['a', 'b']]
        exp_index = CategoricalIndex(['a', 'b'],
                                     categories=index.categories)
        exp = DataFrame({'A': [1, 2], 'B': [4, 5]}, index=exp_index)
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]

        exp_index = CategoricalIndex(['a', 'a', 'b'],
                                     categories=index.categories)
        exp = DataFrame({'A': [1, 1, 2], 'B': [4, 4, 5]}, index=exp_index)
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assert_raises_regex(
                KeyError,
                'a list-indexer must only include values that are '
                'in the categories'):
            df.loc[['a', 'x']]

        # duplicated categories and codes
        index = CategoricalIndex(['a', 'b', 'a'])
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}, index=index)

        # unique slice
        res = df.loc[['a', 'b']]
        exp = DataFrame({'A': [1, 3, 2],
                         'B': [4, 6, 5]},
                        index=CategoricalIndex(['a', 'a', 'b']))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]
        exp = DataFrame(
            {'A': [1, 3, 1, 3, 2],
             'B': [4, 6, 4, 6, 5
                   ]}, index=CategoricalIndex(['a', 'a', 'a', 'a', 'b']))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assert_raises_regex(
                KeyError,
                'a list-indexer must only include values '
                'that are in the categories'):
            df.loc[['a', 'x']]

        # contains unused category
        index = CategoricalIndex(
            ['a', 'b', 'a', 'c'], categories=list('abcde'))
        df = DataFrame({'A': [1, 2, 3, 4], 'B': [5, 6, 7, 8]}, index=index)

        res = df.loc[['a', 'b']]
        exp = DataFrame({'A': [1, 3, 2], 'B': [5, 7, 6]},
                        index=CategoricalIndex(['a', 'a', 'b'],
                                               categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        res = df.loc[['a', 'e']]
        exp = DataFrame({'A': [1, 3, np.nan], 'B': [5, 7, np.nan]},
                        index=CategoricalIndex(['a', 'a', 'e'],
                                               categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]
        exp = DataFrame({'A': [1, 3, 1, 3, 2], 'B': [5, 7, 5, 7, 6]},
                        index=CategoricalIndex(['a', 'a', 'a', 'a', 'b'],
                                               categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assert_raises_regex(
                KeyError,
                'a list-indexer must only include values '
                'that are in the categories'):
            df.loc[['a', 'x']]

    def test_get_indexer_array(self):
        arr = np.array([Timestamp('1999-12-31 00:00:00'),
                        Timestamp('2000-12-31 00:00:00')], dtype=object)
        cats = [Timestamp('1999-12-31 00:00:00'),
                Timestamp('2000-12-31 00:00:00')]
        ci = CategoricalIndex(cats,
                              categories=cats,
                              ordered=False, dtype='category')
        result = ci.get_indexer(arr)
        expected = np.array([0, 1], dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    def test_get_indexer_same_categories_same_order(self):
        ci = CategoricalIndex(['a', 'b'], categories=['a', 'b'])

        result = ci.get_indexer(CategoricalIndex(['b', 'b'],
                                                 categories=['a', 'b']))
        expected = np.array([1, 1], dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    def test_get_indexer_same_categories_different_order(self):
        # https://github.com/pandas-dev/pandas/issues/19551
        ci = CategoricalIndex(['a', 'b'], categories=['a', 'b'])

        result = ci.get_indexer(CategoricalIndex(['b', 'b'],
                                                 categories=['b', 'a']))
        expected = np.array([1, 1], dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    def test_getitem_with_listlike(self):
        # GH 16115
        cats = Categorical([Timestamp('12-31-1999'),
                            Timestamp('12-31-2000')])

        expected = DataFrame([[1, 0], [0, 1]], dtype='uint8',
                             index=[0, 1], columns=cats)
        dummies = pd.get_dummies(cats)
        result = dummies[[c for c in dummies.columns]]
        assert_frame_equal(result, expected)

    def test_setitem_listlike(self):

        # GH 9469
        # properly coerce the input indexers
        np.random.seed(1)
        c = Categorical(np.random.randint(0, 5, size=150000).astype(
            np.int8)).add_categories([-1000])
        indexer = np.array([100000]).astype(np.int64)
        c[indexer] = -1000

        # we are asserting the code result here
        # which maps to the -1000 category
        result = c.codes[np.array([100000]).astype(np.int64)]
        tm.assert_numpy_array_equal(result, np.array([5], dtype='int8'))

    def test_ix_categorical_index(self):
        # GH 12531
        df = DataFrame(np.random.randn(3, 3),
                       index=list('ABC'), columns=list('XYZ'))
        cdf = df.copy()
        cdf.index = CategoricalIndex(df.index)
        cdf.columns = CategoricalIndex(df.columns)

        expect = Series(df.loc['A', :], index=cdf.columns, name='A')
        assert_series_equal(cdf.loc['A', :], expect)

        expect = Series(df.loc[:, 'X'], index=cdf.index, name='X')
        assert_series_equal(cdf.loc[:, 'X'], expect)

        exp_index = CategoricalIndex(list('AB'), categories=['A', 'B', 'C'])
        expect = DataFrame(df.loc[['A', 'B'], :], columns=cdf.columns,
                           index=exp_index)
        assert_frame_equal(cdf.loc[['A', 'B'], :], expect)

        exp_columns = CategoricalIndex(list('XY'),
                                       categories=['X', 'Y', 'Z'])
        expect = DataFrame(df.loc[:, ['X', 'Y']], index=cdf.index,
                           columns=exp_columns)
        assert_frame_equal(cdf.loc[:, ['X', 'Y']], expect)

        # non-unique
        df = DataFrame(np.random.randn(3, 3),
                       index=list('ABA'), columns=list('XYX'))
        cdf = df.copy()
        cdf.index = CategoricalIndex(df.index)
        cdf.columns = CategoricalIndex(df.columns)

        exp_index = CategoricalIndex(list('AA'), categories=['A', 'B'])
        expect = DataFrame(df.loc['A', :], columns=cdf.columns,
                           index=exp_index)
        assert_frame_equal(cdf.loc['A', :], expect)

        exp_columns = CategoricalIndex(list('XX'), categories=['X', 'Y'])
        expect = DataFrame(df.loc[:, 'X'], index=cdf.index,
                           columns=exp_columns)
        assert_frame_equal(cdf.loc[:, 'X'], expect)

        expect = DataFrame(df.loc[['A', 'B'], :], columns=cdf.columns,
                           index=CategoricalIndex(list('AAB')))
        assert_frame_equal(cdf.loc[['A', 'B'], :], expect)

        expect = DataFrame(df.loc[:, ['X', 'Y']], index=cdf.index,
                           columns=CategoricalIndex(list('XXY')))
        assert_frame_equal(cdf.loc[:, ['X', 'Y']], expect)

    def test_read_only_source(self):
        # GH 10043
        rw_array = np.eye(10)
        rw_df = DataFrame(rw_array)

        ro_array = np.eye(10)
        ro_array.setflags(write=False)
        ro_df = DataFrame(ro_array)

        assert_frame_equal(rw_df.iloc[[1, 2, 3]], ro_df.iloc[[1, 2, 3]])
        assert_frame_equal(rw_df.iloc[[1]], ro_df.iloc[[1]])
        assert_series_equal(rw_df.iloc[1], ro_df.iloc[1])
        assert_frame_equal(rw_df.iloc[1:3], ro_df.iloc[1:3])

        assert_frame_equal(rw_df.loc[[1, 2, 3]], ro_df.loc[[1, 2, 3]])
        assert_frame_equal(rw_df.loc[[1]], ro_df.loc[[1]])
        assert_series_equal(rw_df.loc[1], ro_df.loc[1])
        assert_frame_equal(rw_df.loc[1:3], ro_df.loc[1:3])

    def test_reindexing(self):

        # reindexing
        # convert to a regular index
        result = self.df2.reindex(['a', 'b', 'e'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan],
                              'B': Series(list('aaabbe'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3],
                              'B': Series(list('aaabb'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['e'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['d'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['d'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # since we are actually reindexing with a Categorical
        # then return a Categorical
        cats = list('cabe')

        result = self.df2.reindex(Categorical(['a', 'd'], categories=cats))
        expected = DataFrame({'A': [0, 1, 5, np.nan],
                              'B': Series(list('aaad')).astype(
                                  CDT(cats))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(Categorical(['a'], categories=cats))
        expected = DataFrame({'A': [0, 1, 5],
                              'B': Series(list('aaa')).astype(
                                  CDT(cats))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b', 'e'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan],
                              'B': Series(list('aaabbe'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3],
                              'B': Series(list('aaabb'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['e'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # give back the type of categorical that we received
        result = self.df2.reindex(Categorical(
            ['a', 'd'], categories=cats, ordered=True))
        expected = DataFrame(
            {'A': [0, 1, 5, np.nan],
             'B': Series(list('aaad')).astype(
                 CDT(cats, ordered=True))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(Categorical(
            ['a', 'd'], categories=['a', 'd']))
        expected = DataFrame({'A': [0, 1, 5, np.nan],
                              'B': Series(list('aaad')).astype(
                                  CDT(['a', 'd']))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # passed duplicate indexers are not allowed
        pytest.raises(ValueError, lambda: self.df2.reindex(['a', 'a']))

        # args NotImplemented ATM
        pytest.raises(NotImplementedError,
                      lambda: self.df2.reindex(['a'], method='ffill'))
        pytest.raises(NotImplementedError,
                      lambda: self.df2.reindex(['a'], level=1))
        pytest.raises(NotImplementedError,
                      lambda: self.df2.reindex(['a'], limit=2))

    def test_loc_slice(self):
        # slicing
        # not implemented ATM
        # GH9748

        pytest.raises(TypeError, lambda: self.df.loc[1:5])

        # result = df.loc[1:5]
        # expected = df.iloc[[1,2,3,4]]
        # assert_frame_equal(result, expected)

    def test_boolean_selection(self):

        df3 = self.df3
        df4 = self.df4

        result = df3[df3.index == 'a']
        expected = df3.iloc[[]]
        assert_frame_equal(result, expected)

        result = df4[df4.index == 'a']
        expected = df4.iloc[[]]
        assert_frame_equal(result, expected)

        result = df3[df3.index == 1]
        expected = df3.iloc[[0, 1, 3]]
        assert_frame_equal(result, expected)

        result = df4[df4.index == 1]
        expected = df4.iloc[[0, 1, 3]]
        assert_frame_equal(result, expected)

        # since we have an ordered categorical

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=True,
        #         name=u'B')
        result = df3[df3.index < 2]
        expected = df3.iloc[[4]]
        assert_frame_equal(result, expected)

        result = df3[df3.index > 1]
        expected = df3.iloc[[]]
        assert_frame_equal(result, expected)

        # unordered
        # cannot be compared

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=False,
        #         name=u'B')
        pytest.raises(TypeError, lambda: df4[df4.index < 2])
        pytest.raises(TypeError, lambda: df4[df4.index > 1])

    def test_indexing_with_category(self):

        # https://github.com/pandas-dev/pandas/issues/12564
        # consistent result if comparing as Dataframe

        cat = DataFrame({'A': ['foo', 'bar', 'baz']})
        exp = DataFrame({'A': [True, False, False]})

        res = (cat[['A']] == 'foo')
        tm.assert_frame_equal(res, exp)

        cat['A'] = cat['A'].astype('category')

        res = (cat[['A']] == 'foo')
        tm.assert_frame_equal(res, exp)

    def test_map_with_dict_or_series(self):
        orig_values = ['a', 'B', 1, 'a']
        new_values = ['one', 2, 3.0, 'one']
        cur_index = pd.CategoricalIndex(orig_values, name='XXX')
        expected = pd.CategoricalIndex(new_values,
                                       name='XXX', categories=[3.0, 2, 'one'])

        mapper = pd.Series(new_values[:-1], index=orig_values[:-1])
        output = cur_index.map(mapper)
        # Order of categories in output can be different
        tm.assert_index_equal(expected, output)

        mapper = {o: n for o, n in
                  zip(orig_values[:-1], new_values[:-1])}
        output = cur_index.map(mapper)
        # Order of categories in output can be different
        tm.assert_index_equal(expected, output)
