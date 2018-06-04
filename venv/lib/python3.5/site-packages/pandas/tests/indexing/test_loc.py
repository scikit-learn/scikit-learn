""" test label based indexing with loc """

import itertools
import pytest

from warnings import catch_warnings
import numpy as np

import pandas as pd
from pandas.compat import lrange, StringIO
from pandas import Series, DataFrame, Timestamp, date_range, MultiIndex, Index
from pandas.util import testing as tm
from pandas.tests.indexing.common import Base
from pandas.api.types import is_scalar
from pandas.compat import PY2


class TestLoc(Base):

    def test_loc_getitem_dups(self):
        # GH 5678
        # repeated gettitems on a dup index returning a ndarray
        df = DataFrame(
            np.random.random_sample((20, 5)),
            index=['ABCDE' [x % 5] for x in range(20)])
        expected = df.loc['A', 0]
        result = df.loc[:, 0].loc['A']
        tm.assert_series_equal(result, expected)

    def test_loc_getitem_dups2(self):

        # GH4726
        # dup indexing with iloc/loc
        df = DataFrame([[1, 2, 'foo', 'bar', Timestamp('20130101')]],
                       columns=['a', 'a', 'a', 'a', 'a'], index=[1])
        expected = Series([1, 2, 'foo', 'bar', Timestamp('20130101')],
                          index=['a', 'a', 'a', 'a', 'a'], name=1)

        result = df.iloc[0]
        tm.assert_series_equal(result, expected)

        result = df.loc[1]
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_dups(self):

        # GH 6541
        df_orig = DataFrame(
            {'me': list('rttti'),
             'foo': list('aaade'),
             'bar': np.arange(5, dtype='float64') * 1.34 + 2,
             'bar2': np.arange(5, dtype='float64') * -.34 + 2}).set_index('me')

        indexer = tuple(['r', ['bar', 'bar2']])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_series_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

        indexer = tuple(['r', 'bar'])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        assert df.loc[indexer] == 2.0 * df_orig.loc[indexer]

        indexer = tuple(['t', ['bar', 'bar2']])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_frame_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

    def test_loc_setitem_slice(self):
        # GH10503

        # assigning the same type should not change the type
        df1 = DataFrame({'a': [0, 1, 1],
                         'b': Series([100, 200, 300], dtype='uint32')})
        ix = df1['a'] == 1
        newb1 = df1.loc[ix, 'b'] + 1
        df1.loc[ix, 'b'] = newb1
        expected = DataFrame({'a': [0, 1, 1],
                              'b': Series([100, 201, 301], dtype='uint32')})
        tm.assert_frame_equal(df1, expected)

        # assigning a new type should get the inferred type
        df2 = DataFrame({'a': [0, 1, 1], 'b': [100, 200, 300]},
                        dtype='uint64')
        ix = df1['a'] == 1
        newb2 = df2.loc[ix, 'b']
        df1.loc[ix, 'b'] = newb2
        expected = DataFrame({'a': [0, 1, 1], 'b': [100, 200, 300]},
                             dtype='uint64')
        tm.assert_frame_equal(df2, expected)

    def test_loc_getitem_int(self):

        # int label
        self.check_result('int label', 'loc', 2, 'ix', 2,
                          typs=['ints', 'uints'], axes=0)
        self.check_result('int label', 'loc', 3, 'ix', 3,
                          typs=['ints', 'uints'], axes=1)
        self.check_result('int label', 'loc', 4, 'ix', 4,
                          typs=['ints', 'uints'], axes=2)
        self.check_result('int label', 'loc', 2, 'ix', 2,
                          typs=['label'], fails=KeyError)

    def test_loc_getitem_label(self):

        # label
        self.check_result('label', 'loc', 'c', 'ix', 'c', typs=['labels'],
                          axes=0)
        self.check_result('label', 'loc', 'null', 'ix', 'null', typs=['mixed'],
                          axes=0)
        self.check_result('label', 'loc', 8, 'ix', 8, typs=['mixed'], axes=0)
        self.check_result('label', 'loc', Timestamp('20130102'), 'ix', 1,
                          typs=['ts'], axes=0)
        self.check_result('label', 'loc', 'c', 'ix', 'c', typs=['empty'],
                          fails=KeyError)

    def test_loc_getitem_label_out_of_range(self):

        # out of range label
        self.check_result('label range', 'loc', 'f', 'ix', 'f',
                          typs=['ints', 'uints', 'labels', 'mixed', 'ts'],
                          fails=KeyError)
        self.check_result('label range', 'loc', 'f', 'ix', 'f',
                          typs=['floats'], fails=KeyError)
        self.check_result('label range', 'loc', 20, 'ix', 20,
                          typs=['ints', 'uints', 'mixed'], fails=KeyError)
        self.check_result('label range', 'loc', 20, 'ix', 20,
                          typs=['labels'], fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs=['ts'],
                          axes=0, fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs=['floats'],
                          axes=0, fails=KeyError)

    def test_loc_getitem_label_list(self):

        # list of labels
        self.check_result('list lbl', 'loc', [0, 2, 4], 'ix', [0, 2, 4],
                          typs=['ints', 'uints'], axes=0)
        self.check_result('list lbl', 'loc', [3, 6, 9], 'ix', [3, 6, 9],
                          typs=['ints', 'uints'], axes=1)
        self.check_result('list lbl', 'loc', [4, 8, 12], 'ix', [4, 8, 12],
                          typs=['ints', 'uints'], axes=2)
        self.check_result('list lbl', 'loc', ['a', 'b', 'd'], 'ix',
                          ['a', 'b', 'd'], typs=['labels'], axes=0)
        self.check_result('list lbl', 'loc', ['A', 'B', 'C'], 'ix',
                          ['A', 'B', 'C'], typs=['labels'], axes=1)
        self.check_result('list lbl', 'loc', ['Z', 'Y', 'W'], 'ix',
                          ['Z', 'Y', 'W'], typs=['labels'], axes=2)
        self.check_result('list lbl', 'loc', [2, 8, 'null'], 'ix',
                          [2, 8, 'null'], typs=['mixed'], axes=0)
        self.check_result('list lbl', 'loc',
                          [Timestamp('20130102'), Timestamp('20130103')], 'ix',
                          [Timestamp('20130102'), Timestamp('20130103')],
                          typs=['ts'], axes=0)

    @pytest.mark.skipif(PY2, reason=("Catching warnings unreliable with "
                                     "Python 2 (GH #20770)"))
    def test_loc_getitem_label_list_with_missing(self):
        self.check_result('list lbl', 'loc', [0, 1, 2], 'indexer', [0, 1, 2],
                          typs=['empty'], fails=KeyError)
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.check_result('list lbl', 'loc', [0, 2, 10], 'ix', [0, 2, 10],
                              typs=['ints', 'uints', 'floats'],
                              axes=0, fails=KeyError)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.check_result('list lbl', 'loc', [3, 6, 7], 'ix', [3, 6, 7],
                              typs=['ints', 'uints', 'floats'],
                              axes=1, fails=KeyError)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.check_result('list lbl', 'loc', [4, 8, 10], 'ix', [4, 8, 10],
                              typs=['ints', 'uints', 'floats'],
                              axes=2, fails=KeyError)

        # GH 17758 - MultiIndex and missing keys
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.check_result('list lbl', 'loc', [(1, 3), (1, 4), (2, 5)],
                              'ix', [(1, 3), (1, 4), (2, 5)],
                              typs=['multi'],
                              axes=0)

    def test_getitem_label_list_with_missing(self):
        s = Series(range(3), index=['a', 'b', 'c'])

        # consistency
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            s[['a', 'd']]

        s = Series(range(3))
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            s[[0, 3]]

    def test_loc_getitem_label_list_fails(self):
        # fails
        self.check_result('list lbl', 'loc', [20, 30, 40], 'ix', [20, 30, 40],
                          typs=['ints', 'uints'], axes=1, fails=KeyError)
        self.check_result('list lbl', 'loc', [20, 30, 40], 'ix', [20, 30, 40],
                          typs=['ints', 'uints'], axes=2, fails=KeyError)

    def test_loc_getitem_label_array_like(self):
        # array like
        self.check_result('array like', 'loc', Series(index=[0, 2, 4]).index,
                          'ix', [0, 2, 4], typs=['ints', 'uints'], axes=0)
        self.check_result('array like', 'loc', Series(index=[3, 6, 9]).index,
                          'ix', [3, 6, 9], typs=['ints', 'uints'], axes=1)
        self.check_result('array like', 'loc', Series(index=[4, 8, 12]).index,
                          'ix', [4, 8, 12], typs=['ints', 'uints'], axes=2)

    def test_loc_getitem_bool(self):
        # boolean indexers
        b = [True, False, True, False]
        self.check_result('bool', 'loc', b, 'ix', b,
                          typs=['ints', 'uints', 'labels',
                                'mixed', 'ts', 'floats'])
        self.check_result('bool', 'loc', b, 'ix', b, typs=['empty'],
                          fails=KeyError)

    def test_loc_getitem_int_slice(self):

        # ok
        self.check_result('int slice2', 'loc', slice(2, 4), 'ix', [2, 4],
                          typs=['ints', 'uints'], axes=0)
        self.check_result('int slice2', 'loc', slice(3, 6), 'ix', [3, 6],
                          typs=['ints', 'uints'], axes=1)
        self.check_result('int slice2', 'loc', slice(4, 8), 'ix', [4, 8],
                          typs=['ints', 'uints'], axes=2)

        # GH 3053
        # loc should treat integer slices like label slices

        index = MultiIndex.from_tuples([t for t in itertools.product(
            [6, 7, 8], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[6:8, :]
        expected = df
        tm.assert_frame_equal(result, expected)

        index = MultiIndex.from_tuples([t
                                        for t in itertools.product(
                                            [10, 20, 30], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[20:30, :]
        expected = df.iloc[2:]
        tm.assert_frame_equal(result, expected)

        # doc examples
        result = df.loc[10, :]
        expected = df.iloc[0:2]
        expected.index = ['a', 'b']
        tm.assert_frame_equal(result, expected)

        result = df.loc[:, 10]
        # expected = df.ix[:,10] (this fails)
        expected = df[10]
        tm.assert_frame_equal(result, expected)

    def test_loc_to_fail(self):

        # GH3449
        df = DataFrame(np.random.random((3, 3)),
                       index=['a', 'b', 'c'],
                       columns=['e', 'f', 'g'])

        # raise a KeyError?
        pytest.raises(KeyError, df.loc.__getitem__,
                      tuple([[1, 2], [1, 2]]))

        # GH  7496
        # loc should not fallback

        s = Series()
        s.loc[1] = 1
        s.loc['a'] = 2

        pytest.raises(KeyError, lambda: s.loc[-1])
        pytest.raises(KeyError, lambda: s.loc[[-1, -2]])

        pytest.raises(KeyError, lambda: s.loc[['4']])

        s.loc[-1] = 3
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = s.loc[[-1, -2]]
        expected = Series([3, np.nan], index=[-1, -2])
        tm.assert_series_equal(result, expected)

        s['a'] = 2
        pytest.raises(KeyError, lambda: s.loc[[-2]])

        del s['a']

        def f():
            s.loc[[-2]] = 0

        pytest.raises(KeyError, f)

        # inconsistency between .loc[values] and .loc[values,:]
        # GH 7999
        df = DataFrame([['a'], ['b']], index=[1, 2], columns=['value'])

        def f():
            df.loc[[3], :]

        pytest.raises(KeyError, f)

        def f():
            df.loc[[3]]

        pytest.raises(KeyError, f)

    def test_loc_getitem_list_with_fail(self):
        # 15747
        # should KeyError if *any* missing labels

        s = Series([1, 2, 3])

        s.loc[[2]]

        with pytest.raises(KeyError):
            s.loc[[3]]

        # a non-match and a match
        with tm.assert_produces_warning(FutureWarning):
            expected = s.loc[[2, 3]]
        result = s.reindex([2, 3])
        tm.assert_series_equal(result, expected)

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)
        self.check_result('lab slice', 'loc', slice(1, 3),
                          'ix', slice(1, 3),
                          typs=['labels', 'mixed', 'empty', 'ts', 'floats'],
                          fails=TypeError)

        # real label slices
        self.check_result('lab slice', 'loc', slice('a', 'c'),
                          'ix', slice('a', 'c'), typs=['labels'], axes=0)
        self.check_result('lab slice', 'loc', slice('A', 'C'),
                          'ix', slice('A', 'C'), typs=['labels'], axes=1)
        self.check_result('lab slice', 'loc', slice('W', 'Z'),
                          'ix', slice('W', 'Z'), typs=['labels'], axes=2)

        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=0)
        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=1, fails=TypeError)
        self.check_result('ts  slice', 'loc', slice('20130102', '20130104'),
                          'ix', slice('20130102', '20130104'),
                          typs=['ts'], axes=2, fails=TypeError)

        # GH 14316
        self.check_result('ts slice rev', 'loc', slice('20130104', '20130102'),
                          'indexer', [0, 1, 2], typs=['ts_rev'], axes=0)

        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=0, fails=TypeError)
        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=1, fails=KeyError)
        self.check_result('mixed slice', 'loc', slice(2, 8), 'ix', slice(2, 8),
                          typs=['mixed'], axes=2, fails=KeyError)

        self.check_result('mixed slice', 'loc', slice(2, 4, 2), 'ix', slice(
            2, 4, 2), typs=['mixed'], axes=0, fails=TypeError)

    def test_loc_index(self):
        # gh-17131
        # a boolean index should index like a boolean numpy array

        df = DataFrame(
            np.random.random(size=(5, 10)),
            index=["alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1"])

        mask = df.index.map(lambda x: "alpha" in x)
        expected = df.loc[np.array(mask)]

        result = df.loc[mask]
        tm.assert_frame_equal(result, expected)

        result = df.loc[mask.values]
        tm.assert_frame_equal(result, expected)

    def test_loc_general(self):

        df = DataFrame(
            np.random.rand(4, 4), columns=['A', 'B', 'C', 'D'],
            index=['A', 'B', 'C', 'D'])

        # want this to work
        result = df.loc[:, "A":"B"].iloc[0:2, :]
        assert (result.columns == ['A', 'B']).all()
        assert (result.index == ['A', 'B']).all()

        # mixed type
        result = DataFrame({'a': [Timestamp('20130101')], 'b': [1]}).iloc[0]
        expected = Series([Timestamp('20130101'), 1], index=['a', 'b'], name=0)
        tm.assert_series_equal(result, expected)
        assert result.dtype == object

    def test_loc_setitem_consistency(self):
        # GH 6149
        # coerce similarly for setitem and loc when rows have a null-slice
        expected = DataFrame({'date': Series(0, index=range(5),
                                             dtype=np.int64),
                              'val': Series(range(5), dtype=np.int64)})

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(
                            range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 0
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = np.array(0, dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = np.array([0, 0, 0, 0, 0], dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        expected = DataFrame({'date': Series('foo', index=range(5)),
                              'val': Series(range(5), dtype=np.int64)})
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 'foo'
        tm.assert_frame_equal(df, expected)

        expected = DataFrame({'date': Series(1.0, index=range(5)),
                              'val': Series(range(5), dtype=np.int64)})
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(range(5), dtype=np.int64)})
        df.loc[:, 'date'] = 1.0
        tm.assert_frame_equal(df, expected)

        # GH 15494
        # setting on frame with single row
        df = DataFrame({'date': Series([Timestamp('20180101')])})
        df.loc[:, 'date'] = 'string'
        expected = DataFrame({'date': Series(['string'])})
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_empty(self):
        # empty (essentially noops)
        expected = DataFrame(columns=['x', 'y'])
        expected['x'] = expected['x'].astype(np.int64)
        df = DataFrame(columns=['x', 'y'])
        df.loc[:, 'x'] = 1
        tm.assert_frame_equal(df, expected)

        df = DataFrame(columns=['x', 'y'])
        df['x'] = 1
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_slice_column_len(self):
        # .loc[:,column] setting with slice == len of the column
        # GH10408
        data = """Level_0,,,Respondent,Respondent,Respondent,OtherCat,OtherCat
Level_1,,,Something,StartDate,EndDate,Yes/No,SomethingElse
Region,Site,RespondentID,,,,,
Region_1,Site_1,3987227376,A,5/25/2015 10:59,5/25/2015 11:22,Yes,
Region_1,Site_1,3980680971,A,5/21/2015 9:40,5/21/2015 9:52,Yes,Yes
Region_1,Site_2,3977723249,A,5/20/2015 8:27,5/20/2015 8:41,Yes,
Region_1,Site_2,3977723089,A,5/20/2015 8:33,5/20/2015 9:09,Yes,No"""

        df = pd.read_csv(StringIO(data), header=[0, 1], index_col=[0, 1, 2])
        df.loc[:, ('Respondent', 'StartDate')] = pd.to_datetime(df.loc[:, (
            'Respondent', 'StartDate')])
        df.loc[:, ('Respondent', 'EndDate')] = pd.to_datetime(df.loc[:, (
            'Respondent', 'EndDate')])
        df.loc[:, ('Respondent', 'Duration')] = df.loc[:, (
            'Respondent', 'EndDate')] - df.loc[:, ('Respondent', 'StartDate')]

        df.loc[:, ('Respondent', 'Duration')] = df.loc[:, (
            'Respondent', 'Duration')].astype('timedelta64[s]')
        expected = Series([1380, 720, 840, 2160.], index=df.index,
                          name=('Respondent', 'Duration'))
        tm.assert_series_equal(df[('Respondent', 'Duration')], expected)

    def test_loc_setitem_frame(self):
        df = self.frame_labels

        result = df.iloc[0, 0]

        df.loc['a', 'A'] = 1
        result = df.loc['a', 'A']
        assert result == 1

        result = df.iloc[0, 0]
        assert result == 1

        df.loc[:, 'B':'D'] = 0
        expected = df.loc[:, 'B':'D']
        result = df.iloc[:, 1:]
        tm.assert_frame_equal(result, expected)

        # GH 6254
        # setting issue
        df = DataFrame(index=[3, 5, 4], columns=['A'])
        df.loc[[4, 3, 5], 'A'] = np.array([1, 2, 3], dtype='int64')
        expected = DataFrame(dict(A=Series(
            [1, 2, 3], index=[4, 3, 5]))).reindex(index=[3, 5, 4])
        tm.assert_frame_equal(df, expected)

        # GH 6252
        # setting with an empty frame
        keys1 = ['@' + str(i) for i in range(5)]
        val1 = np.arange(5, dtype='int64')

        keys2 = ['@' + str(i) for i in range(4)]
        val2 = np.arange(4, dtype='int64')

        index = list(set(keys1).union(keys2))
        df = DataFrame(index=index)
        df['A'] = np.nan
        df.loc[keys1, 'A'] = val1

        df['B'] = np.nan
        df.loc[keys2, 'B'] = val2

        expected = DataFrame(dict(A=Series(val1, index=keys1), B=Series(
            val2, index=keys2))).reindex(index=index)
        tm.assert_frame_equal(df, expected)

        # GH 8669
        # invalid coercion of nan -> int
        df = DataFrame({'A': [1, 2, 3], 'B': np.nan})
        df.loc[df.B > df.A, 'B'] = df.A
        expected = DataFrame({'A': [1, 2, 3], 'B': np.nan})
        tm.assert_frame_equal(df, expected)

        # GH 6546
        # setting with mixed labels
        df = DataFrame({1: [1, 2], 2: [3, 4], 'a': ['a', 'b']})

        result = df.loc[0, [1, 2]]
        expected = Series([1, 3], index=[1, 2], dtype=object, name=0)
        tm.assert_series_equal(result, expected)

        expected = DataFrame({1: [5, 2], 2: [6, 4], 'a': ['a', 'b']})
        df.loc[0, [1, 2]] = [5, 6]
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_frame_multiples(self):
        # multiple setting
        df = DataFrame({'A': ['foo', 'bar', 'baz'],
                        'B': Series(
                            range(3), dtype=np.int64)})
        rhs = df.loc[1:2]
        rhs.index = df.index[0:2]
        df.loc[0:1] = rhs
        expected = DataFrame({'A': ['bar', 'baz', 'baz'],
                              'B': Series(
                                  [1, 2, 2], dtype=np.int64)})
        tm.assert_frame_equal(df, expected)

        # multiple setting with frame on rhs (with M8)
        df = DataFrame({'date': date_range('2000-01-01', '2000-01-5'),
                        'val': Series(
                            range(5), dtype=np.int64)})
        expected = DataFrame({'date': [Timestamp('20000101'), Timestamp(
            '20000102'), Timestamp('20000101'), Timestamp('20000102'),
            Timestamp('20000103')],
            'val': Series(
            [0, 1, 0, 1, 2], dtype=np.int64)})
        rhs = df.loc[0:2]
        rhs.index = df.index[2:5]
        df.loc[2:4] = rhs
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        'indexer', [['A'], slice(None, 'A', None), np.array(['A'])])
    @pytest.mark.parametrize(
        'value', [['Z'], np.array(['Z'])])
    def test_loc_setitem_with_scalar_index(self, indexer, value):
        # GH #19474
        # assigning like "df.loc[0, ['A']] = ['Z']" should be evaluated
        # elementwisely, not using "setter('A', ['Z'])".

        df = pd.DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        df.loc[0, indexer] = value
        result = df.loc[0, 'A']

        assert is_scalar(result) and result == 'Z'

    def test_loc_coerceion(self):

        # 12411
        df = DataFrame({'date': [Timestamp('20130101').tz_localize('UTC'),
                                 pd.NaT]})
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 12045
        import datetime
        df = DataFrame({'date': [datetime.datetime(2012, 1, 1),
                                 datetime.datetime(1012, 1, 2)]})
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 11594
        df = DataFrame({'text': ['some words'] + [None] * 9})
        expected = df.dtypes

        result = df.iloc[0:2]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[3:]
        tm.assert_series_equal(result.dtypes, expected)

    def test_loc_non_unique(self):
        # GH3659
        # non-unique indexer with loc slice
        # https://groups.google.com/forum/?fromgroups#!topic/pydata/zTm2No0crYs

        # these are going to raise because the we are non monotonic
        df = DataFrame({'A': [1, 2, 3, 4, 5, 6],
                        'B': [3, 4, 5, 6, 7, 8]}, index=[0, 1, 0, 1, 2, 3])
        pytest.raises(KeyError, df.loc.__getitem__,
                      tuple([slice(1, None)]))
        pytest.raises(KeyError, df.loc.__getitem__,
                      tuple([slice(0, None)]))
        pytest.raises(KeyError, df.loc.__getitem__, tuple([slice(1, 2)]))

        # monotonic are ok
        df = DataFrame({'A': [1, 2, 3, 4, 5, 6],
                        'B': [3, 4, 5, 6, 7, 8]},
                       index=[0, 1, 0, 1, 2, 3]).sort_index(axis=0)
        result = df.loc[1:]
        expected = DataFrame({'A': [2, 4, 5, 6], 'B': [4, 6, 7, 8]},
                             index=[1, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

        result = df.loc[0:]
        tm.assert_frame_equal(result, df)

        result = df.loc[1:2]
        expected = DataFrame({'A': [2, 4, 5], 'B': [4, 6, 7]},
                             index=[1, 1, 2])
        tm.assert_frame_equal(result, expected)

    def test_loc_non_unique_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list('ABCDEFG')

        def gen_test(l, l2):
            return pd.concat([
                DataFrame(np.random.randn(l, len(columns)),
                          index=lrange(l), columns=columns),
                DataFrame(np.ones((l2, len(columns))),
                          index=[0] * l2, columns=columns)])

        def gen_expected(df, mask):
            l = len(mask)
            return pd.concat([df.take([0]),
                              DataFrame(np.ones((l, len(columns))),
                                        index=[0] * l,
                                        columns=columns),
                              df.take(mask[1:])])

        df = gen_test(900, 100)
        assert not df.index.is_unique

        mask = np.arange(100)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

        df = gen_test(900000, 100000)
        assert not df.index.is_unique

        mask = np.arange(100000)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

    def test_loc_name(self):
        # GH 3880
        df = DataFrame([[1, 1], [1, 1]])
        df.index.name = 'index_name'
        result = df.iloc[[0, 1]].index.name
        assert result == 'index_name'

        with catch_warnings(record=True):
            result = df.ix[[0, 1]].index.name
        assert result == 'index_name'

        result = df.loc[[0, 1]].index.name
        assert result == 'index_name'

    def test_loc_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        tm.assert_frame_equal(df.loc[:, []], df.iloc[:, :0],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.loc[[], :], df.iloc[:0, :],
                              check_index_type=True, check_column_type=True)
        # horizontal empty
        tm.assert_frame_equal(df.loc[[]], df.iloc[:0, :],
                              check_index_type=True,
                              check_column_type=True)

    def test_identity_slice_returns_new_object(self):
        # GH13873
        original_df = DataFrame({'a': [1, 2, 3]})
        sliced_df = original_df.loc[:]
        assert sliced_df is not original_df
        assert original_df[:] is not original_df

        # should be a shallow copy
        original_df['a'] = [4, 4, 4]
        assert (sliced_df['a'] == 4).all()

        # These should not return copies
        assert original_df is original_df.loc[:, :]
        df = DataFrame(np.random.randn(10, 4))
        assert df[0] is df.loc[:, 0]

        # Same tests for Series
        original_series = Series([1, 2, 3, 4, 5, 6])
        sliced_series = original_series.loc[:]
        assert sliced_series is not original_series
        assert original_series[:] is not original_series

        original_series[:3] = [7, 8, 9]
        assert all(sliced_series[:3] == [7, 8, 9])

    @pytest.mark.parametrize(
        'indexer_type_1',
        (list, tuple, set, slice, np.ndarray, Series, Index))
    @pytest.mark.parametrize(
        'indexer_type_2',
        (list, tuple, set, slice, np.ndarray, Series, Index))
    def test_loc_getitem_nested_indexer(self, indexer_type_1, indexer_type_2):
        # GH #19686
        # .loc should work with nested indexers which can be
        # any list-like objects (see `pandas.api.types.is_list_like`) or slices

        def convert_nested_indexer(indexer_type, keys):
            if indexer_type == np.ndarray:
                return np.array(keys)
            if indexer_type == slice:
                return slice(*keys)
            return indexer_type(keys)

        a = [10, 20, 30]
        b = [1, 2, 3]
        index = pd.MultiIndex.from_product([a, b])
        df = pd.DataFrame(
            np.arange(len(index), dtype='int64'),
            index=index, columns=['Data'])

        keys = ([10, 20], [2, 3])
        types = (indexer_type_1, indexer_type_2)

        # check indexers with all the combinations of nested objects
        # of all the valid types
        indexer = tuple(
            convert_nested_indexer(indexer_type, k)
            for indexer_type, k in zip(types, keys))

        result = df.loc[indexer, 'Data']
        expected = pd.Series(
            [1, 2, 4, 5], name='Data',
            index=pd.MultiIndex.from_product(keys))

        tm.assert_series_equal(result, expected)

    def test_loc_uint64(self):
        # GH20722
        # Test whether loc accept uint64 max value as index.
        s = pd.Series([1, 2],
                      index=[np.iinfo('uint64').max - 1,
                             np.iinfo('uint64').max])

        result = s.loc[np.iinfo('uint64').max - 1]
        expected = s.iloc[0]
        assert result == expected

        result = s.loc[[np.iinfo('uint64').max - 1]]
        expected = s.iloc[[0]]
        tm.assert_series_equal(result, expected)

        result = s.loc[[np.iinfo('uint64').max - 1,
                       np.iinfo('uint64').max]]
        tm.assert_series_equal(result, s)
