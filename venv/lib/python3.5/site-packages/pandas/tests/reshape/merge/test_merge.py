# pylint: disable=E1103

import random
import re
from collections import OrderedDict
from datetime import date, datetime

import numpy as np
import pytest
from numpy import nan
from numpy.random import randn

import pandas as pd
import pandas.util.testing as tm
from pandas import (Categorical, CategoricalIndex, DataFrame, DatetimeIndex,
                    Float64Index, Index, Int64Index, MultiIndex, RangeIndex,
                    Series, UInt64Index)
from pandas.api.types import CategoricalDtype as CDT
from pandas.compat import lrange, lzip
from pandas.core.dtypes.common import is_categorical_dtype, is_object_dtype
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.core.reshape.concat import concat
from pandas.core.reshape.merge import MergeError, merge
from pandas.util.testing import assert_frame_equal, assert_series_equal

N = 50
NGROUPS = 8


def get_test_data(ngroups=NGROUPS, n=N):
    unique_groups = lrange(ngroups)
    arr = np.asarray(np.tile(unique_groups, n // ngroups))

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)])

    random.shuffle(arr)
    return arr


class TestMerge(object):

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

        self.left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                               'v1': np.random.randn(7)})
        self.right = DataFrame({'v2': np.random.randn(4)},
                               index=['d', 'b', 'c', 'a'])

    def test_merge_inner_join_empty(self):
        # GH 15328
        df_empty = pd.DataFrame()
        df_a = pd.DataFrame({'a': [1, 2]}, index=[0, 1], dtype='int64')
        result = pd.merge(df_empty, df_a, left_index=True, right_index=True)
        expected = pd.DataFrame({'a': []}, index=[], dtype='int64')
        assert_frame_equal(result, expected)

    def test_merge_common(self):
        joined = merge(self.df, self.df2)
        exp = merge(self.df, self.df2, on=['key1', 'key2'])
        tm.assert_frame_equal(joined, exp)

    def test_merge_index_as_on_arg(self):
        # GH14355

        left = self.df.set_index('key1')
        right = self.df2.set_index('key1')
        result = merge(left, right, on='key1')
        expected = merge(self.df, self.df2, on='key1').set_index('key1')
        assert_frame_equal(result, expected)

    def test_merge_index_singlekey_right_vs_left(self):
        left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                          'v1': np.random.randn(7)})
        right = DataFrame({'v2': np.random.randn(4)},
                          index=['d', 'b', 'c', 'a'])

        merged1 = merge(left, right, left_on='key',
                        right_index=True, how='left', sort=False)
        merged2 = merge(right, left, right_on='key',
                        left_index=True, how='right', sort=False)
        assert_frame_equal(merged1, merged2.loc[:, merged1.columns])

        merged1 = merge(left, right, left_on='key',
                        right_index=True, how='left', sort=True)
        merged2 = merge(right, left, right_on='key',
                        left_index=True, how='right', sort=True)
        assert_frame_equal(merged1, merged2.loc[:, merged1.columns])

    def test_merge_index_singlekey_inner(self):
        left = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'e', 'a'],
                          'v1': np.random.randn(7)})
        right = DataFrame({'v2': np.random.randn(4)},
                          index=['d', 'b', 'c', 'a'])

        # inner join
        result = merge(left, right, left_on='key', right_index=True,
                       how='inner')
        expected = left.join(right, on='key').loc[result.index]
        assert_frame_equal(result, expected)

        result = merge(right, left, right_on='key', left_index=True,
                       how='inner')
        expected = left.join(right, on='key').loc[result.index]
        assert_frame_equal(result, expected.loc[:, result.columns])

    def test_merge_misspecified(self):
        pytest.raises(ValueError, merge, self.left, self.right,
                      left_index=True)
        pytest.raises(ValueError, merge, self.left, self.right,
                      right_index=True)

        pytest.raises(ValueError, merge, self.left, self.left,
                      left_on='key', on='key')

        pytest.raises(ValueError, merge, self.df, self.df2,
                      left_on=['key1'], right_on=['key1', 'key2'])

    def test_index_and_on_parameters_confusion(self):
        pytest.raises(ValueError, merge, self.df, self.df2, how='left',
                      left_index=False, right_index=['key1', 'key2'])
        pytest.raises(ValueError, merge, self.df, self.df2, how='left',
                      left_index=['key1', 'key2'], right_index=False)
        pytest.raises(ValueError, merge, self.df, self.df2, how='left',
                      left_index=['key1', 'key2'],
                      right_index=['key1', 'key2'])

    def test_merge_overlap(self):
        merged = merge(self.left, self.left, on='key')
        exp_len = (self.left['key'].value_counts() ** 2).sum()
        assert len(merged) == exp_len
        assert 'v1_x' in merged
        assert 'v1_y' in merged

    def test_merge_different_column_key_names(self):
        left = DataFrame({'lkey': ['foo', 'bar', 'baz', 'foo'],
                          'value': [1, 2, 3, 4]})
        right = DataFrame({'rkey': ['foo', 'bar', 'qux', 'foo'],
                           'value': [5, 6, 7, 8]})

        merged = left.merge(right, left_on='lkey', right_on='rkey',
                            how='outer', sort=True)

        exp = pd.Series(['bar', 'baz', 'foo', 'foo', 'foo', 'foo', np.nan],
                        name='lkey')
        tm.assert_series_equal(merged['lkey'], exp)

        exp = pd.Series(['bar', np.nan, 'foo', 'foo', 'foo', 'foo', 'qux'],
                        name='rkey')
        tm.assert_series_equal(merged['rkey'], exp)

        exp = pd.Series([2, 3, 1, 1, 4, 4, np.nan], name='value_x')
        tm.assert_series_equal(merged['value_x'], exp)

        exp = pd.Series([6, np.nan, 5, 8, 5, 8, 7], name='value_y')
        tm.assert_series_equal(merged['value_y'], exp)

    def test_merge_copy(self):
        left = DataFrame({'a': 0, 'b': 1}, index=lrange(10))
        right = DataFrame({'c': 'foo', 'd': 'bar'}, index=lrange(10))

        merged = merge(left, right, left_index=True,
                       right_index=True, copy=True)

        merged['a'] = 6
        assert (left['a'] == 0).all()

        merged['d'] = 'peekaboo'
        assert (right['d'] == 'bar').all()

    def test_merge_nocopy(self):
        left = DataFrame({'a': 0, 'b': 1}, index=lrange(10))
        right = DataFrame({'c': 'foo', 'd': 'bar'}, index=lrange(10))

        merged = merge(left, right, left_index=True,
                       right_index=True, copy=False)

        merged['a'] = 6
        assert (left['a'] == 6).all()

        merged['d'] = 'peekaboo'
        assert (right['d'] == 'peekaboo').all()

    def test_intelligently_handle_join_key(self):
        # #733, be a bit more 1337 about not returning unconsolidated DataFrame

        left = DataFrame({'key': [1, 1, 2, 2, 3],
                          'value': lrange(5)}, columns=['value', 'key'])
        right = DataFrame({'key': [1, 1, 2, 3, 4, 5],
                           'rvalue': lrange(6)})

        joined = merge(left, right, on='key', how='outer')
        expected = DataFrame({'key': [1, 1, 1, 1, 2, 2, 3, 4, 5],
                              'value': np.array([0, 0, 1, 1, 2, 3, 4,
                                                 np.nan, np.nan]),
                              'rvalue': [0, 1, 0, 1, 2, 2, 3, 4, 5]},
                             columns=['value', 'key', 'rvalue'])
        assert_frame_equal(joined, expected)

    def test_merge_join_key_dtype_cast(self):
        # #8596

        df1 = DataFrame({'key': [1], 'v1': [10]})
        df2 = DataFrame({'key': [2], 'v1': [20]})
        df = merge(df1, df2, how='outer')
        assert df['key'].dtype == 'int64'

        df1 = DataFrame({'key': [True], 'v1': [1]})
        df2 = DataFrame({'key': [False], 'v1': [0]})
        df = merge(df1, df2, how='outer')

        # GH13169
        # this really should be bool
        assert df['key'].dtype == 'object'

        df1 = DataFrame({'val': [1]})
        df2 = DataFrame({'val': [2]})
        lkey = np.array([1])
        rkey = np.array([2])
        df = merge(df1, df2, left_on=lkey, right_on=rkey, how='outer')
        assert df['key_0'].dtype == 'int64'

    def test_handle_join_key_pass_array(self):
        left = DataFrame({'key': [1, 1, 2, 2, 3],
                          'value': lrange(5)}, columns=['value', 'key'])
        right = DataFrame({'rvalue': lrange(6)})
        key = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on='key', right_on=key, how='outer')
        merged2 = merge(right, left, left_on=key, right_on='key', how='outer')

        assert_series_equal(merged['key'], merged2['key'])
        assert merged['key'].notna().all()
        assert merged2['key'].notna().all()

        left = DataFrame({'value': lrange(5)}, columns=['value'])
        right = DataFrame({'rvalue': lrange(6)})
        lkey = np.array([1, 1, 2, 2, 3])
        rkey = np.array([1, 1, 2, 3, 4, 5])

        merged = merge(left, right, left_on=lkey, right_on=rkey, how='outer')
        tm.assert_series_equal(merged['key_0'], Series([1, 1, 1, 1, 2,
                                                        2, 3, 4, 5],
                                                       name='key_0'))

        left = DataFrame({'value': lrange(3)})
        right = DataFrame({'rvalue': lrange(6)})

        key = np.array([0, 1, 1, 2, 2, 3], dtype=np.int64)
        merged = merge(left, right, left_index=True, right_on=key, how='outer')
        tm.assert_series_equal(merged['key_0'], Series(key, name='key_0'))

    def test_no_overlap_more_informative_error(self):
        dt = datetime.now()
        df1 = DataFrame({'x': ['a']}, index=[dt])

        df2 = DataFrame({'y': ['b', 'c']}, index=[dt, dt])
        pytest.raises(MergeError, merge, df1, df2)

        msg = ('No common columns to perform merge on. '
               'Merge options: left_on={lon}, right_on={ron}, '
               'left_index={lidx}, right_index={ridx}'
               .format(lon=None, ron=None, lidx=False, ridx=False))

        with tm.assert_raises_regex(MergeError, msg):
            merge(df1, df2)

    def test_merge_non_unique_indexes(self):

        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        dt4 = datetime(2012, 5, 4)

        df1 = DataFrame({'x': ['a']}, index=[dt])
        df2 = DataFrame({'y': ['b', 'c']}, index=[dt, dt])
        _check_merge(df1, df2)

        # Not monotonic
        df1 = DataFrame({'x': ['a', 'b', 'q']}, index=[dt2, dt, dt4])
        df2 = DataFrame({'y': ['c', 'd', 'e', 'f', 'g', 'h']},
                        index=[dt3, dt3, dt2, dt2, dt, dt])
        _check_merge(df1, df2)

        df1 = DataFrame({'x': ['a', 'b']}, index=[dt, dt])
        df2 = DataFrame({'y': ['c', 'd']}, index=[dt, dt])
        _check_merge(df1, df2)

    def test_merge_non_unique_index_many_to_many(self):
        dt = datetime(2012, 5, 1)
        dt2 = datetime(2012, 5, 2)
        dt3 = datetime(2012, 5, 3)
        df1 = DataFrame({'x': ['a', 'b', 'c', 'd']},
                        index=[dt2, dt2, dt, dt])
        df2 = DataFrame({'y': ['e', 'f', 'g', ' h', 'i']},
                        index=[dt2, dt2, dt3, dt, dt])
        _check_merge(df1, df2)

    def test_left_merge_empty_dataframe(self):
        left = DataFrame({'key': [1], 'value': [2]})
        right = DataFrame({'key': []})

        result = merge(left, right, on='key', how='left')
        assert_frame_equal(result, left)

        result = merge(right, left, on='key', how='right')
        assert_frame_equal(result, left)

    @pytest.mark.parametrize('kwarg',
                             [dict(left_index=True, right_index=True),
                              dict(left_index=True, right_on='x'),
                              dict(left_on='a', right_index=True),
                              dict(left_on='a', right_on='x')])
    def test_merge_left_empty_right_empty(self, join_type, kwarg):
        # GH 10824
        left = pd.DataFrame([], columns=['a', 'b', 'c'])
        right = pd.DataFrame([], columns=['x', 'y', 'z'])

        exp_in = pd.DataFrame([], columns=['a', 'b', 'c', 'x', 'y', 'z'],
                              index=pd.Index([], dtype=object),
                              dtype=object)

        result = pd.merge(left, right, how=join_type, **kwarg)
        tm.assert_frame_equal(result, exp_in)

    def test_merge_left_empty_right_notempty(self):
        # GH 10824
        left = pd.DataFrame([], columns=['a', 'b', 'c'])
        right = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                             columns=['x', 'y', 'z'])

        exp_out = pd.DataFrame({'a': np.array([np.nan] * 3, dtype=object),
                                'b': np.array([np.nan] * 3, dtype=object),
                                'c': np.array([np.nan] * 3, dtype=object),
                                'x': [1, 4, 7],
                                'y': [2, 5, 8],
                                'z': [3, 6, 9]},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp)

        for kwarg in [dict(left_index=True, right_index=True),
                      dict(left_index=True, right_on='x')]:
            check1(exp_in, kwarg)
            check2(exp_out, kwarg)

        kwarg = dict(left_on='a', right_index=True)
        check1(exp_in, kwarg)
        exp_out['a'] = [0, 1, 2]
        check2(exp_out, kwarg)

        kwarg = dict(left_on='a', right_on='x')
        check1(exp_in, kwarg)
        exp_out['a'] = np.array([np.nan] * 3, dtype=object)
        check2(exp_out, kwarg)

    def test_merge_left_notempty_right_empty(self):
        # GH 10824
        left = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            columns=['a', 'b', 'c'])
        right = pd.DataFrame([], columns=['x', 'y', 'z'])

        exp_out = pd.DataFrame({'a': [1, 4, 7],
                                'b': [2, 5, 8],
                                'c': [3, 6, 9],
                                'x': np.array([np.nan] * 3, dtype=object),
                                'y': np.array([np.nan] * 3, dtype=object),
                                'z': np.array([np.nan] * 3, dtype=object)},
                               columns=['a', 'b', 'c', 'x', 'y', 'z'])
        exp_in = exp_out[0:0]  # make empty DataFrame keeping dtype
        # result will have object dtype
        exp_in.index = exp_in.index.astype(object)

        def check1(exp, kwarg):
            result = pd.merge(left, right, how='inner', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='right', **kwarg)
            tm.assert_frame_equal(result, exp)

        def check2(exp, kwarg):
            result = pd.merge(left, right, how='left', **kwarg)
            tm.assert_frame_equal(result, exp)
            result = pd.merge(left, right, how='outer', **kwarg)
            tm.assert_frame_equal(result, exp)

            for kwarg in [dict(left_index=True, right_index=True),
                          dict(left_index=True, right_on='x'),
                          dict(left_on='a', right_index=True),
                          dict(left_on='a', right_on='x')]:
                check1(exp_in, kwarg)
                check2(exp_out, kwarg)

    def test_merge_nosort(self):
        # #2098, anything to do?

        from datetime import datetime

        d = {"var1": np.random.randint(0, 10, size=10),
             "var2": np.random.randint(0, 10, size=10),
             "var3": [datetime(2012, 1, 12),
                      datetime(2011, 2, 4),
                      datetime(2010, 2, 3),
                      datetime(2012, 1, 12),
                      datetime(2011, 2, 4),
                      datetime(2012, 4, 3),
                      datetime(2012, 3, 4),
                      datetime(2008, 5, 1),
                      datetime(2010, 2, 3),
                      datetime(2012, 2, 3)]}
        df = DataFrame.from_dict(d)
        var3 = df.var3.unique()
        var3.sort()
        new = DataFrame.from_dict({"var3": var3,
                                   "var8": np.random.random(7)})

        result = df.merge(new, on="var3", sort=False)
        exp = merge(df, new, on='var3', sort=False)
        assert_frame_equal(result, exp)

        assert (df.var3.unique() == result.var3.unique()).all()

    def test_merge_nan_right(self):
        df1 = DataFrame({"i1": [0, 1], "i2": [0, 1]})
        df2 = DataFrame({"i1": [0], "i3": [0]})
        result = df1.join(df2, on="i1", rsuffix="_")
        expected = (DataFrame({'i1': {0: 0.0, 1: 1}, 'i2': {0: 0, 1: 1},
                               'i1_': {0: 0, 1: np.nan},
                               'i3': {0: 0.0, 1: np.nan},
                               None: {0: 0, 1: 0}})
                    .set_index(None)
                    .reset_index()[['i1', 'i2', 'i1_', 'i3']])
        assert_frame_equal(result, expected, check_dtype=False)

        df1 = DataFrame({"i1": [0, 1], "i2": [0.5, 1.5]})
        df2 = DataFrame({"i1": [0], "i3": [0.7]})
        result = df1.join(df2, rsuffix="_", on='i1')
        expected = (DataFrame({'i1': {0: 0, 1: 1}, 'i1_': {0: 0.0, 1: nan},
                               'i2': {0: 0.5, 1: 1.5},
                               'i3': {0: 0.69999999999999996,
                                      1: nan}})
                    [['i1', 'i2', 'i1_', 'i3']])
        assert_frame_equal(result, expected)

    def test_merge_type(self):
        class NotADataFrame(DataFrame):

            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.df)
        result = nad.merge(self.df2, on='key1')

        assert isinstance(result, NotADataFrame)

    def test_join_append_timedeltas(self):

        import datetime as dt
        from pandas import NaT

        # timedelta64 issues with join/merge
        # GH 5695

        d = {'d': dt.datetime(2013, 11, 5, 5, 56), 't': dt.timedelta(0, 22500)}
        df = DataFrame(columns=list('dt'))
        df = df.append(d, ignore_index=True)
        result = df.append(d, ignore_index=True)
        expected = DataFrame({'d': [dt.datetime(2013, 11, 5, 5, 56),
                                    dt.datetime(2013, 11, 5, 5, 56)],
                              't': [dt.timedelta(0, 22500),
                                    dt.timedelta(0, 22500)]})
        assert_frame_equal(result, expected)

        td = np.timedelta64(300000000)
        lhs = DataFrame(Series([td, td], index=["A", "B"]))
        rhs = DataFrame(Series([td], index=["A"]))

        result = lhs.join(rhs, rsuffix='r', how="left")
        expected = DataFrame({'0': Series([td, td], index=list('AB')),
                              '0r': Series([td, NaT], index=list('AB'))})
        assert_frame_equal(result, expected)

    def test_other_datetime_unit(self):
        # GH 13389
        df1 = pd.DataFrame({'entity_id': [101, 102]})
        s = pd.Series([None, None], index=[101, 102], name='days')

        for dtype in ['datetime64[D]', 'datetime64[h]', 'datetime64[m]',
                      'datetime64[s]', 'datetime64[ms]', 'datetime64[us]',
                      'datetime64[ns]']:

            df2 = s.astype(dtype).to_frame('days')
            # coerces to datetime64[ns], thus sholuld not be affected
            assert df2['days'].dtype == 'datetime64[ns]'

            result = df1.merge(df2, left_on='entity_id', right_index=True)

            exp = pd.DataFrame({'entity_id': [101, 102],
                                'days': np.array(['nat', 'nat'],
                                                 dtype='datetime64[ns]')},
                               columns=['entity_id', 'days'])
            tm.assert_frame_equal(result, exp)

    @pytest.mark.parametrize("unit", ['D', 'h', 'm', 's', 'ms', 'us', 'ns'])
    def test_other_timedelta_unit(self, unit):
        # GH 13389
        df1 = pd.DataFrame({'entity_id': [101, 102]})
        s = pd.Series([None, None], index=[101, 102], name='days')

        dtype = "m8[{}]".format(unit)
        df2 = s.astype(dtype).to_frame('days')
        assert df2['days'].dtype == 'm8[ns]'

        result = df1.merge(df2, left_on='entity_id', right_index=True)

        exp = pd.DataFrame({'entity_id': [101, 102],
                            'days': np.array(['nat', 'nat'],
                                             dtype=dtype)},
                           columns=['entity_id', 'days'])
        tm.assert_frame_equal(result, exp)

    def test_overlapping_columns_error_message(self):
        df = DataFrame({'key': [1, 2, 3],
                        'v1': [4, 5, 6],
                        'v2': [7, 8, 9]})
        df2 = DataFrame({'key': [1, 2, 3],
                         'v1': [4, 5, 6],
                         'v2': [7, 8, 9]})

        df.columns = ['key', 'foo', 'foo']
        df2.columns = ['key', 'bar', 'bar']
        expected = DataFrame({'key': [1, 2, 3],
                              'v1': [4, 5, 6],
                              'v2': [7, 8, 9],
                              'v3': [4, 5, 6],
                              'v4': [7, 8, 9]})
        expected.columns = ['key', 'foo', 'foo', 'bar', 'bar']
        assert_frame_equal(merge(df, df2), expected)

        # #2649, #10639
        df2.columns = ['key1', 'foo', 'foo']
        pytest.raises(ValueError, merge, df, df2)

    def test_merge_on_datetime64tz(self):

        # GH11405
        left = pd.DataFrame({'key': pd.date_range('20151010', periods=2,
                                                  tz='US/Eastern'),
                             'value': [1, 2]})
        right = pd.DataFrame({'key': pd.date_range('20151011', periods=3,
                                                   tz='US/Eastern'),
                              'value': [1, 2, 3]})

        expected = DataFrame({'key': pd.date_range('20151010', periods=4,
                                                   tz='US/Eastern'),
                              'value_x': [1, 2, np.nan, np.nan],
                              'value_y': [np.nan, 1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)

        left = pd.DataFrame({'key': [1, 2],
                             'value': pd.date_range('20151010', periods=2,
                                                    tz='US/Eastern')})
        right = pd.DataFrame({'key': [2, 3],
                              'value': pd.date_range('20151011', periods=2,
                                                     tz='US/Eastern')})
        expected = DataFrame({
            'key': [1, 2, 3],
            'value_x': list(pd.date_range('20151010', periods=2,
                                          tz='US/Eastern')) + [pd.NaT],
            'value_y': [pd.NaT] + list(pd.date_range('20151011', periods=2,
                                                     tz='US/Eastern'))})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)
        assert result['value_x'].dtype == 'datetime64[ns, US/Eastern]'
        assert result['value_y'].dtype == 'datetime64[ns, US/Eastern]'

    def test_merge_non_unique_period_index(self):
        # GH #16871
        index = pd.period_range('2016-01-01', periods=16, freq='M')
        df = DataFrame([i for i in range(len(index))],
                       index=index, columns=['pnum'])
        df2 = concat([df, df])
        result = df.merge(df2, left_index=True, right_index=True, how='inner')
        expected = DataFrame(
            np.tile(np.arange(16, dtype=np.int64).repeat(2).reshape(-1, 1), 2),
            columns=['pnum_x', 'pnum_y'], index=df2.sort_index().index)
        tm.assert_frame_equal(result, expected)

    def test_merge_on_periods(self):
        left = pd.DataFrame({'key': pd.period_range('20151010', periods=2,
                                                    freq='D'),
                             'value': [1, 2]})
        right = pd.DataFrame({'key': pd.period_range('20151011', periods=3,
                                                     freq='D'),
                              'value': [1, 2, 3]})

        expected = DataFrame({'key': pd.period_range('20151010', periods=4,
                                                     freq='D'),
                              'value_x': [1, 2, np.nan, np.nan],
                              'value_y': [np.nan, 1, 2, 3]})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)

        left = pd.DataFrame({'key': [1, 2],
                             'value': pd.period_range('20151010', periods=2,
                                                      freq='D')})
        right = pd.DataFrame({'key': [2, 3],
                              'value': pd.period_range('20151011', periods=2,
                                                       freq='D')})

        exp_x = pd.period_range('20151010', periods=2, freq='D')
        exp_y = pd.period_range('20151011', periods=2, freq='D')
        expected = DataFrame({'key': [1, 2, 3],
                              'value_x': list(exp_x) + [pd.NaT],
                              'value_y': [pd.NaT] + list(exp_y)})
        result = pd.merge(left, right, on='key', how='outer')
        assert_frame_equal(result, expected)
        assert result['value_x'].dtype == 'object'
        assert result['value_y'].dtype == 'object'

    def test_indicator(self):
        # PR #10054. xref #7412 and closes #8790.
        df1 = DataFrame({'col1': [0, 1], 'col_conflict': [1, 2],
                         'col_left': ['a', 'b']})
        df1_copy = df1.copy()

        df2 = DataFrame({'col1': [1, 2, 3, 4, 5],
                         'col_conflict': [1, 2, 3, 4, 5],
                         'col_right': [2, 2, 2, 2, 2]})
        df2_copy = df2.copy()

        df_result = DataFrame({
            'col1': [0, 1, 2, 3, 4, 5],
            'col_conflict_x': [1, 2, np.nan, np.nan, np.nan, np.nan],
            'col_left': ['a', 'b', np.nan, np.nan, np.nan, np.nan],
            'col_conflict_y': [np.nan, 1, 2, 3, 4, 5],
            'col_right': [np.nan, 2, 2, 2, 2, 2]})
        df_result['_merge'] = Categorical(
            ['left_only', 'both', 'right_only',
             'right_only', 'right_only', 'right_only'],
            categories=['left_only', 'right_only', 'both'])

        df_result = df_result[['col1', 'col_conflict_x', 'col_left',
                               'col_conflict_y', 'col_right', '_merge']]

        test = merge(df1, df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)
        test = df1.merge(df2, on='col1', how='outer', indicator=True)
        assert_frame_equal(test, df_result)

        # No side effects
        assert_frame_equal(df1, df1_copy)
        assert_frame_equal(df2, df2_copy)

        # Check with custom name
        df_result_custom_name = df_result
        df_result_custom_name = df_result_custom_name.rename(
            columns={'_merge': 'custom_name'})

        test_custom_name = merge(
            df1, df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)
        test_custom_name = df1.merge(
            df2, on='col1', how='outer', indicator='custom_name')
        assert_frame_equal(test_custom_name, df_result_custom_name)

        # Check only accepts strings and booleans
        with pytest.raises(ValueError):
            merge(df1, df2, on='col1', how='outer', indicator=5)
        with pytest.raises(ValueError):
            df1.merge(df2, on='col1', how='outer', indicator=5)

        # Check result integrity

        test2 = merge(df1, df2, on='col1', how='left', indicator=True)
        assert (test2._merge != 'right_only').all()
        test2 = df1.merge(df2, on='col1', how='left', indicator=True)
        assert (test2._merge != 'right_only').all()

        test3 = merge(df1, df2, on='col1', how='right', indicator=True)
        assert (test3._merge != 'left_only').all()
        test3 = df1.merge(df2, on='col1', how='right', indicator=True)
        assert (test3._merge != 'left_only').all()

        test4 = merge(df1, df2, on='col1', how='inner', indicator=True)
        assert (test4._merge == 'both').all()
        test4 = df1.merge(df2, on='col1', how='inner', indicator=True)
        assert (test4._merge == 'both').all()

        # Check if working name in df
        for i in ['_right_indicator', '_left_indicator', '_merge']:
            df_badcolumn = DataFrame({'col1': [1, 2], i: [2, 2]})

            with pytest.raises(ValueError):
                merge(df1, df_badcolumn, on='col1',
                      how='outer', indicator=True)
            with pytest.raises(ValueError):
                df1.merge(df_badcolumn, on='col1', how='outer', indicator=True)

        # Check for name conflict with custom name
        df_badcolumn = DataFrame(
            {'col1': [1, 2], 'custom_column_name': [2, 2]})

        with pytest.raises(ValueError):
            merge(df1, df_badcolumn, on='col1', how='outer',
                  indicator='custom_column_name')
        with pytest.raises(ValueError):
            df1.merge(df_badcolumn, on='col1', how='outer',
                      indicator='custom_column_name')

        # Merge on multiple columns
        df3 = DataFrame({'col1': [0, 1], 'col2': ['a', 'b']})

        df4 = DataFrame({'col1': [1, 1, 3], 'col2': ['b', 'x', 'y']})

        hand_coded_result = DataFrame({'col1': [0, 1, 1, 3],
                                       'col2': ['a', 'b', 'x', 'y']})
        hand_coded_result['_merge'] = Categorical(
            ['left_only', 'both', 'right_only', 'right_only'],
            categories=['left_only', 'right_only', 'both'])

        test5 = merge(df3, df4, on=['col1', 'col2'],
                      how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)
        test5 = df3.merge(df4, on=['col1', 'col2'],
                          how='outer', indicator=True)
        assert_frame_equal(test5, hand_coded_result)

    def test_validation(self):
        left = DataFrame({'a': ['a', 'b', 'c', 'd'],
                          'b': ['cat', 'dog', 'weasel', 'horse']},
                         index=range(4))

        right = DataFrame({'a': ['a', 'b', 'c', 'd', 'e'],
                           'c': ['meow', 'bark', 'um... weasel noise?',
                                 'nay', 'chirp']},
                          index=range(5))

        # Make sure no side effects.
        left_copy = left.copy()
        right_copy = right.copy()

        result = merge(left, right, left_index=True, right_index=True,
                       validate='1:1')
        assert_frame_equal(left, left_copy)
        assert_frame_equal(right, right_copy)

        # make sure merge still correct
        expected = DataFrame({'a_x': ['a', 'b', 'c', 'd'],
                              'b': ['cat', 'dog', 'weasel', 'horse'],
                              'a_y': ['a', 'b', 'c', 'd'],
                              'c': ['meow', 'bark', 'um... weasel noise?',
                                    'nay']},
                             index=range(4),
                             columns=['a_x', 'b', 'a_y', 'c'])

        result = merge(left, right, left_index=True, right_index=True,
                       validate='one_to_one')
        assert_frame_equal(result, expected)

        expected_2 = DataFrame({'a': ['a', 'b', 'c', 'd'],
                                'b': ['cat', 'dog', 'weasel', 'horse'],
                                'c': ['meow', 'bark', 'um... weasel noise?',
                                      'nay']},
                               index=range(4))

        result = merge(left, right, on='a', validate='1:1')
        assert_frame_equal(left, left_copy)
        assert_frame_equal(right, right_copy)
        assert_frame_equal(result, expected_2)

        result = merge(left, right, on='a', validate='one_to_one')
        assert_frame_equal(result, expected_2)

        # One index, one column
        expected_3 = DataFrame({'b': ['cat', 'dog', 'weasel', 'horse'],
                                'a': ['a', 'b', 'c', 'd'],
                                'c': ['meow', 'bark', 'um... weasel noise?',
                                      'nay']},
                               columns=['b', 'a', 'c'],
                               index=range(4))

        left_index_reset = left.set_index('a')
        result = merge(left_index_reset, right, left_index=True,
                       right_on='a', validate='one_to_one')
        assert_frame_equal(result, expected_3)

        # Dups on right
        right_w_dups = right.append(pd.DataFrame({'a': ['e'], 'c': ['moo']},
                                                 index=[4]))
        merge(left, right_w_dups, left_index=True, right_index=True,
              validate='one_to_many')

        with pytest.raises(MergeError):
            merge(left, right_w_dups, left_index=True, right_index=True,
                  validate='one_to_one')

        with pytest.raises(MergeError):
            merge(left, right_w_dups, on='a', validate='one_to_one')

        # Dups on left
        left_w_dups = left.append(pd.DataFrame({'a': ['a'], 'c': ['cow']},
                                               index=[3]), sort=True)
        merge(left_w_dups, right, left_index=True, right_index=True,
              validate='many_to_one')

        with pytest.raises(MergeError):
            merge(left_w_dups, right, left_index=True, right_index=True,
                  validate='one_to_one')

        with pytest.raises(MergeError):
            merge(left_w_dups, right, on='a', validate='one_to_one')

        # Dups on both
        merge(left_w_dups, right_w_dups, on='a', validate='many_to_many')

        with pytest.raises(MergeError):
            merge(left_w_dups, right_w_dups, left_index=True,
                  right_index=True, validate='many_to_one')

        with pytest.raises(MergeError):
            merge(left_w_dups, right_w_dups, on='a',
                  validate='one_to_many')

        # Check invalid arguments
        with pytest.raises(ValueError):
            merge(left, right, on='a', validate='jibberish')

        # Two column merge, dups in both, but jointly no dups.
        left = DataFrame({'a': ['a', 'a', 'b', 'b'],
                          'b': [0, 1, 0, 1],
                          'c': ['cat', 'dog', 'weasel', 'horse']},
                         index=range(4))

        right = DataFrame({'a': ['a', 'a', 'b'],
                           'b': [0, 1, 0],
                           'd': ['meow', 'bark', 'um... weasel noise?']},
                          index=range(3))

        expected_multi = DataFrame({'a': ['a', 'a', 'b'],
                                    'b': [0, 1, 0],
                                    'c': ['cat', 'dog', 'weasel'],
                                    'd': ['meow', 'bark',
                                          'um... weasel noise?']},
                                   index=range(3))

        with pytest.raises(MergeError):
            merge(left, right, on='a', validate='1:1')

        result = merge(left, right, on=['a', 'b'], validate='1:1')
        assert_frame_equal(result, expected_multi)

    def test_merge_two_empty_df_no_division_error(self):
        # GH17776, PR #17846
        a = pd.DataFrame({'a': [], 'b': [], 'c': []})
        with np.errstate(divide='raise'):
            merge(a, a, on=('a', 'b'))


def _check_merge(x, y):
    for how in ['inner', 'left', 'outer']:
        result = x.join(y, how=how)

        expected = merge(x.reset_index(), y.reset_index(), how=how,
                         sort=True)
        expected = expected.set_index('index')

        # TODO check_names on merge?
        assert_frame_equal(result, expected, check_names=False)


class TestMergeMulti(object):

    def setup_method(self, method):
        self.index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                        ['one', 'two', 'three']],
                                labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                        [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                                names=['first', 'second'])
        self.to_join = DataFrame(np.random.randn(10, 3), index=self.index,
                                 columns=['j_one', 'j_two', 'j_three'])

        # a little relevant example with NAs
        key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
                'qux', 'snap']
        key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
                'three', 'one']

        data = np.random.randn(len(key1))
        self.data = DataFrame({'key1': key1, 'key2': key2,
                               'data': data})

    def test_merge_on_multikey(self):
        joined = self.data.join(self.to_join, on=['key1', 'key2'])

        join_key = Index(lzip(self.data['key1'], self.data['key2']))
        indexer = self.to_join.index.get_indexer(join_key)
        ex_values = self.to_join.values.take(indexer, axis=0)
        ex_values[indexer == -1] = np.nan
        expected = self.data.join(DataFrame(ex_values,
                                            columns=self.to_join.columns))

        # TODO: columns aren't in the same order yet
        assert_frame_equal(joined, expected.loc[:, joined.columns])

        left = self.data.join(self.to_join, on=['key1', 'key2'], sort=True)
        right = expected.loc[:, joined.columns].sort_values(['key1', 'key2'],
                                                            kind='mergesort')
        assert_frame_equal(left, right)

    def test_left_join_multi_index(self):
        icols = ['1st', '2nd', '3rd']

        def bind_cols(df):
            iord = lambda a: 0 if a != a else ord(a)
            f = lambda ts: ts.map(iord) - ord('a')
            return (f(df['1st']) + f(df['3rd']) * 1e2 +
                    df['2nd'].fillna(0) * 1e4)

        def run_asserts(left, right):
            for sort in [False, True]:
                res = left.join(right, on=icols, how='left', sort=sort)

                assert len(left) < len(res) + 1
                assert not res['4th'].isna().any()
                assert not res['5th'].isna().any()

                tm.assert_series_equal(
                    res['4th'], - res['5th'], check_names=False)
                result = bind_cols(res.iloc[:, :-2])
                tm.assert_series_equal(res['4th'], result, check_names=False)
                assert result.name is None

                if sort:
                    tm.assert_frame_equal(
                        res, res.sort_values(icols, kind='mergesort'))

                out = merge(left, right.reset_index(), on=icols,
                            sort=sort, how='left')

                res.index = np.arange(len(res))
                tm.assert_frame_equal(out, res)

        lc = list(map(chr, np.arange(ord('a'), ord('z') + 1)))
        left = DataFrame(np.random.choice(lc, (5000, 2)),
                         columns=['1st', '3rd'])
        left.insert(1, '2nd', np.random.randint(0, 1000, len(left)))

        i = np.random.permutation(len(left))
        right = left.iloc[i].copy()

        left['4th'] = bind_cols(left)
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)

        # inject some nulls
        left.loc[1::23, '1st'] = np.nan
        left.loc[2::37, '2nd'] = np.nan
        left.loc[3::43, '3rd'] = np.nan
        left['4th'] = bind_cols(left)

        i = np.random.permutation(len(left))
        right = left.iloc[i, :-1]
        right['5th'] = - bind_cols(right)
        right.set_index(icols, inplace=True)

        run_asserts(left, right)

    def test_merge_right_vs_left(self):
        # compare left vs right merge with multikey
        for sort in [False, True]:
            merged1 = self.data.merge(self.to_join, left_on=['key1', 'key2'],
                                      right_index=True, how='left', sort=sort)

            merged2 = self.to_join.merge(self.data, right_on=['key1', 'key2'],
                                         left_index=True, how='right',
                                         sort=sort)

            merged2 = merged2.loc[:, merged1.columns]
            assert_frame_equal(merged1, merged2)

    def test_compress_group_combinations(self):

        # ~ 40000000 possible unique groups
        key1 = tm.rands_array(10, 10000)
        key1 = np.tile(key1, 2)
        key2 = key1[::-1]

        df = DataFrame({'key1': key1, 'key2': key2,
                        'value1': np.random.randn(20000)})

        df2 = DataFrame({'key1': key1[::2], 'key2': key2[::2],
                         'value2': np.random.randn(10000)})

        # just to hit the label compression code path
        merge(df, df2, how='outer')

    def test_left_join_index_preserve_order(self):

        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'v': np.array(np.arange(24), dtype=np.int64)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(
            result.sort_values(['k1', 'k2'], kind='mergesort'),
            left.join(right, on=['k1', 'k2'], sort=True))

        # test join with multi dtypes blocks
        left = DataFrame({'k1': [0, 1, 2] * 8,
                          'k2': ['foo', 'bar'] * 12,
                          'k3': np.array([0, 1, 2] * 8, dtype=np.float32),
                          'v': np.array(np.arange(24), dtype=np.int32)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': [5, 7]}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()
        expected['v2'] = np.nan
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(
            result.sort_values(['k1', 'k2'], kind='mergesort'),
            left.join(right, on=['k1', 'k2'], sort=True))

        # do a right join for an extra test
        joined = merge(right, left, left_index=True,
                       right_on=['k1', 'k2'], how='right')
        tm.assert_frame_equal(joined.loc[:, expected.columns], expected)

    def test_left_join_index_multi_match_multiindex(self):
        left = DataFrame([
            ['X', 'Y', 'C', 'a'],
            ['W', 'Y', 'C', 'e'],
            ['V', 'Q', 'A', 'h'],
            ['V', 'R', 'D', 'i'],
            ['X', 'Y', 'D', 'b'],
            ['X', 'Y', 'A', 'c'],
            ['W', 'Q', 'B', 'f'],
            ['W', 'R', 'C', 'g'],
            ['V', 'Y', 'C', 'j'],
            ['X', 'Y', 'B', 'd']],
            columns=['cola', 'colb', 'colc', 'tag'],
            index=[3, 2, 0, 1, 7, 6, 4, 5, 9, 8])

        right = DataFrame([
            ['W', 'R', 'C', 0],
            ['W', 'Q', 'B', 3],
            ['W', 'Q', 'B', 8],
            ['X', 'Y', 'A', 1],
            ['X', 'Y', 'A', 4],
            ['X', 'Y', 'B', 5],
            ['X', 'Y', 'C', 6],
            ['X', 'Y', 'C', 9],
            ['X', 'Q', 'C', -6],
            ['X', 'R', 'C', -9],
            ['V', 'Y', 'C', 7],
            ['V', 'R', 'D', 2],
            ['V', 'R', 'D', -1],
            ['V', 'Q', 'A', -3]],
            columns=['col1', 'col2', 'col3', 'val'])

        right.set_index(['col1', 'col2', 'col3'], inplace=True)
        result = left.join(right, on=['cola', 'colb', 'colc'], how='left')

        expected = DataFrame([
            ['X', 'Y', 'C', 'a', 6],
            ['X', 'Y', 'C', 'a', 9],
            ['W', 'Y', 'C', 'e', nan],
            ['V', 'Q', 'A', 'h', -3],
            ['V', 'R', 'D', 'i', 2],
            ['V', 'R', 'D', 'i', -1],
            ['X', 'Y', 'D', 'b', nan],
            ['X', 'Y', 'A', 'c', 1],
            ['X', 'Y', 'A', 'c', 4],
            ['W', 'Q', 'B', 'f', 3],
            ['W', 'Q', 'B', 'f', 8],
            ['W', 'R', 'C', 'g', 0],
            ['V', 'Y', 'C', 'j', 7],
            ['X', 'Y', 'B', 'd', 5]],
            columns=['cola', 'colb', 'colc', 'tag', 'val'],
            index=[3, 3, 2, 0, 1, 1, 7, 6, 6, 4, 4, 5, 9, 8])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=['cola', 'colb', 'colc'],
                           how='left', sort=True)

        tm.assert_frame_equal(
            result,
            expected.sort_values(['cola', 'colb', 'colc'], kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        right.reset_index(inplace=True)
        right.columns = left.columns[:3].tolist() + right.columns[-1:].tolist()
        result = merge(left, right, how='left', on=left.columns[:-1].tolist())
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_left_join_index_multi_match(self):
        left = DataFrame([
            ['c', 0],
            ['b', 1],
            ['a', 2],
            ['b', 3]],
            columns=['tag', 'val'],
            index=[2, 0, 1, 3])

        right = DataFrame([
            ['a', 'v'],
            ['c', 'w'],
            ['c', 'x'],
            ['d', 'y'],
            ['a', 'z'],
            ['c', 'r'],
            ['e', 'q'],
            ['c', 's']],
            columns=['tag', 'char'])

        right.set_index('tag', inplace=True)
        result = left.join(right, on='tag', how='left')

        expected = DataFrame([
            ['c', 0, 'w'],
            ['c', 0, 'x'],
            ['c', 0, 'r'],
            ['c', 0, 's'],
            ['b', 1, nan],
            ['a', 2, 'v'],
            ['a', 2, 'z'],
            ['b', 3, nan]],
            columns=['tag', 'val', 'char'],
            index=[2, 2, 2, 2, 0, 1, 1, 3])

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on='tag', how='left', sort=True)
        tm.assert_frame_equal(
            result, expected.sort_values('tag', kind='mergesort'))

        # GH7331 - maintain left frame order in left merge
        result = merge(left, right.reset_index(), how='left', on='tag')
        expected.index = np.arange(len(expected))
        tm.assert_frame_equal(result, expected)

    def test_left_merge_na_buglet(self):
        left = DataFrame({'id': list('abcde'), 'v1': randn(5),
                          'v2': randn(5), 'dummy': list('abcde'),
                          'v3': randn(5)},
                         columns=['id', 'v1', 'v2', 'dummy', 'v3'])
        right = DataFrame({'id': ['a', 'b', np.nan, np.nan, np.nan],
                           'sv3': [1.234, 5.678, np.nan, np.nan, np.nan]})

        merged = merge(left, right, on='id', how='left')

        rdf = right.drop(['id'], axis=1)
        expected = left.join(rdf)
        tm.assert_frame_equal(merged, expected)

    def test_merge_na_keys(self):
        data = [[1950, "A", 1.5],
                [1950, "B", 1.5],
                [1955, "B", 1.5],
                [1960, "B", np.nan],
                [1970, "B", 4.],
                [1950, "C", 4.],
                [1960, "C", np.nan],
                [1965, "C", 3.],
                [1970, "C", 4.]]

        frame = DataFrame(data, columns=["year", "panel", "data"])

        other_data = [[1960, 'A', np.nan],
                      [1970, 'A', np.nan],
                      [1955, 'A', np.nan],
                      [1965, 'A', np.nan],
                      [1965, 'B', np.nan],
                      [1955, 'C', np.nan]]
        other = DataFrame(other_data, columns=['year', 'panel', 'data'])

        result = frame.merge(other, how='outer')

        expected = frame.fillna(-999).merge(other.fillna(-999), how='outer')
        expected = expected.replace(-999, np.nan)

        tm.assert_frame_equal(result, expected)

    def test_join_multi_levels(self):

        # GH 3662
        # merge multi-levels
        household = (
            DataFrame(
                dict(household_id=[1, 2, 3],
                     male=[0, 1, 0],
                     wealth=[196087.3, 316478.7, 294750]),
                columns=['household_id', 'male', 'wealth'])
            .set_index('household_id'))
        portfolio = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000289783", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     name=["ABN Amro", "Robeco", "Royal Dutch Shell",
                           "Royal Dutch Shell",
                           "AAB Eastern Europe Equity Fund",
                           "Postbank BioTech Fonds", np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'name', 'share'])
            .set_index(['household_id', 'asset_id']))
        result = household.join(portfolio, how='inner')
        expected = (
            DataFrame(
                dict(male=[0, 1, 1, 0, 0, 0],
                     wealth=[196087.3, 316478.7, 316478.7,
                             294750.0, 294750.0, 294750.0],
                     name=['ABN Amro', 'Robeco', 'Royal Dutch Shell',
                           'Royal Dutch Shell',
                           'AAB Eastern Europe Equity Fund',
                           'Postbank BioTech Fonds'],
                     share=[1.00, 0.40, 0.60, 0.15, 0.60, 0.25],
                     household_id=[1, 2, 2, 3, 3, 3],
                     asset_id=['nl0000301109', 'nl0000289783', 'gb00b03mlx29',
                               'gb00b03mlx29', 'lu0197800237',
                               'nl0000289965']))
            .set_index(['household_id', 'asset_id'])
            .reindex(columns=['male', 'wealth', 'name', 'share']))
        assert_frame_equal(result, expected)

        assert_frame_equal(result, expected)

        # equivalency
        result2 = (merge(household.reset_index(), portfolio.reset_index(),
                         on=['household_id'], how='inner')
                   .set_index(['household_id', 'asset_id']))
        assert_frame_equal(result2, expected)

        result = household.join(portfolio, how='outer')
        expected = (concat([
            expected,
            (DataFrame(
                dict(share=[1.00]),
                index=MultiIndex.from_tuples(
                    [(4, np.nan)],
                    names=['household_id', 'asset_id'])))
        ], axis=0, sort=True).reindex(columns=expected.columns))
        assert_frame_equal(result, expected)

        # invalid cases
        household.index.name = 'foo'

        def f():
            household.join(portfolio, how='inner')

        pytest.raises(ValueError, f)

        portfolio2 = portfolio.copy()
        portfolio2.index.set_names(['household_id', 'foo'])

        def f():
            portfolio2.join(portfolio, how='inner')

        pytest.raises(ValueError, f)

    def test_join_multi_levels2(self):

        # some more advanced merges
        # GH6360
        household = (
            DataFrame(
                dict(household_id=[1, 2, 2, 3, 3, 3, 4],
                     asset_id=["nl0000301109", "nl0000301109", "gb00b03mlx29",
                               "gb00b03mlx29", "lu0197800237", "nl0000289965",
                               np.nan],
                     share=[1.0, 0.4, 0.6, 0.15, 0.6, 0.25, 1.0]),
                columns=['household_id', 'asset_id', 'share'])
            .set_index(['household_id', 'asset_id']))

        log_return = DataFrame(dict(
            asset_id=["gb00b03mlx29", "gb00b03mlx29",
                      "gb00b03mlx29", "lu0197800237", "lu0197800237"],
            t=[233, 234, 235, 180, 181],
            log_return=[.09604978, -.06524096, .03532373, .03025441, .036997]
        )).set_index(["asset_id", "t"])

        expected = (
            DataFrame(dict(
                household_id=[2, 2, 2, 3, 3, 3, 3, 3],
                asset_id=["gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237"],
                t=[233, 234, 235, 233, 234, 235, 180, 181],
                share=[0.6, 0.6, 0.6, 0.15, 0.15, 0.15, 0.6, 0.6],
                log_return=[.09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997]
            ))
            .set_index(["household_id", "asset_id", "t"])
            .reindex(columns=['share', 'log_return']))

        def f():
            household.join(log_return, how='inner')

        pytest.raises(NotImplementedError, f)

        # this is the equivalency
        result = (merge(household.reset_index(), log_return.reset_index(),
                        on=['asset_id'], how='inner')
                  .set_index(['household_id', 'asset_id', 't']))
        assert_frame_equal(result, expected)

        expected = (
            DataFrame(dict(
                household_id=[1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4],
                asset_id=["nl0000301109", "nl0000289783", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29",
                          "gb00b03mlx29", "gb00b03mlx29", "gb00b03mlx29",
                          "lu0197800237", "lu0197800237",
                          "nl0000289965", None],
                t=[None, None, 233, 234, 235, 233, 234,
                   235, 180, 181, None, None],
                share=[1.0, 0.4, 0.6, 0.6, 0.6, 0.15,
                       0.15, 0.15, 0.6, 0.6, 0.25, 1.0],
                log_return=[None, None, .09604978, -.06524096, .03532373,
                            .09604978, -.06524096, .03532373,
                            .03025441, .036997, None, None]
            ))
            .set_index(["household_id", "asset_id", "t"]))

        def f():
            household.join(log_return, how='outer')

        pytest.raises(NotImplementedError, f)

    @pytest.mark.parametrize("klass", [None, np.asarray, Series, Index])
    def test_merge_datetime_index(self, klass):
        # see gh-19038
        df = DataFrame([1, 2, 3],
                       ["2016-01-01", "2017-01-01", "2018-01-01"],
                       columns=["a"])
        df.index = pd.to_datetime(df.index)
        on_vector = df.index.year

        if klass is not None:
            on_vector = klass(on_vector)

        expected = DataFrame(
            OrderedDict([
                ("a", [1, 2, 3]),
                ("key_1", [2016, 2017, 2018]),
            ])
        )

        result = df.merge(df, on=["a", on_vector], how="inner")
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            OrderedDict([
                ("key_0", [2016, 2017, 2018]),
                ("a_x", [1, 2, 3]),
                ("a_y", [1, 2, 3]),
            ])
        )

        result = df.merge(df, on=[df.index.year], how="inner")
        tm.assert_frame_equal(result, expected)


class TestMergeDtypes(object):

    @pytest.mark.parametrize('right_vals', [
        ['foo', 'bar'],
        Series(['foo', 'bar']).astype('category'),
        [1, 2],
        [1.0, 2.0],
        Series([1, 2], dtype='uint64'),
        Series([1, 2], dtype='int32')
    ])
    def test_different(self, right_vals):

        left = DataFrame({'A': ['foo', 'bar'],
                          'B': Series(['foo', 'bar']).astype('category'),
                          'C': [1, 2],
                          'D': [1.0, 2.0],
                          'E': Series([1, 2], dtype='uint64'),
                          'F': Series([1, 2], dtype='int32')})
        right = DataFrame({'A': right_vals})

        # GH 9780
        # We allow merging on object and categorical cols and cast
        # categorical cols to object
        if (is_categorical_dtype(right['A'].dtype) or
                is_object_dtype(right['A'].dtype)):
            result = pd.merge(left, right, on='A')
            assert is_object_dtype(result.A.dtype)

        # GH 9780
        # We raise for merging on object col and int/float col and
        # merging on categorical col and int/float col
        else:
            msg = ("You are trying to merge on "
                   "{lk_dtype} and {rk_dtype} columns. "
                   "If you wish to proceed you should use "
                   "pd.concat".format(lk_dtype=left['A'].dtype,
                                      rk_dtype=right['A'].dtype))
            with tm.assert_raises_regex(ValueError, msg):
                pd.merge(left, right, on='A')

    @pytest.mark.parametrize('d1', [np.int64, np.int32,
                                    np.int16, np.int8, np.uint8])
    @pytest.mark.parametrize('d2', [np.int64, np.float64,
                                    np.float32, np.float16])
    def test_join_multi_dtypes(self, d1, d2):

        dtype1 = np.dtype(d1)
        dtype2 = np.dtype(d2)

        left = DataFrame({'k1': np.array([0, 1, 2] * 8, dtype=dtype1),
                          'k2': ['foo', 'bar'] * 12,
                          'v': np.array(np.arange(24), dtype=np.int64)})

        index = MultiIndex.from_tuples([(2, 'bar'), (1, 'foo')])
        right = DataFrame({'v2': np.array([5, 7], dtype=dtype2)}, index=index)

        result = left.join(right, on=['k1', 'k2'])

        expected = left.copy()

        if dtype2.kind == 'i':
            dtype2 = np.dtype('float64')
        expected['v2'] = np.array(np.nan, dtype=dtype2)
        expected.loc[(expected.k1 == 2) & (expected.k2 == 'bar'), 'v2'] = 5
        expected.loc[(expected.k1 == 1) & (expected.k2 == 'foo'), 'v2'] = 7

        tm.assert_frame_equal(result, expected)

        result = left.join(right, on=['k1', 'k2'], sort=True)
        expected.sort_values(['k1', 'k2'], kind='mergesort', inplace=True)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('int_vals, float_vals, exp_vals', [
        ([1, 2, 3], [1.0, 2.0, 3.0], {'X': [1, 2, 3], 'Y': [1.0, 2.0, 3.0]}),
        ([1, 2, 3], [1.0, 3.0], {'X': [1, 3], 'Y': [1.0, 3.0]}),
        ([1, 2], [1.0, 2.0, 3.0], {'X': [1, 2], 'Y': [1.0, 2.0]}),
    ])
    def test_merge_on_ints_floats(self, int_vals, float_vals, exp_vals):
        # GH 16572
        # Check that float column is not cast to object if
        # merging on float and int columns
        A = DataFrame({'X': int_vals})
        B = DataFrame({'Y': float_vals})
        expected = DataFrame(exp_vals)

        result = A.merge(B, left_on='X', right_on='Y')
        assert_frame_equal(result, expected)

        result = B.merge(A, left_on='Y', right_on='X')
        assert_frame_equal(result, expected[['Y', 'X']])

    def test_merge_on_ints_floats_warning(self):
        # GH 16572
        # merge will produce a warning when merging on int and
        # float columns where the float values are not exactly
        # equal to their int representation
        A = DataFrame({'X': [1, 2, 3]})
        B = DataFrame({'Y': [1.1, 2.5, 3.0]})
        expected = DataFrame({'X': [3], 'Y': [3.0]})

        with tm.assert_produces_warning(UserWarning):
            result = A.merge(B, left_on='X', right_on='Y')
            assert_frame_equal(result, expected)

        with tm.assert_produces_warning(UserWarning):
            result = B.merge(A, left_on='Y', right_on='X')
            assert_frame_equal(result, expected[['Y', 'X']])

        # test no warning if float has NaNs
        B = DataFrame({'Y': [np.nan, np.nan, 3.0]})

        with tm.assert_produces_warning(None):
            result = B.merge(A, left_on='Y', right_on='X')
            assert_frame_equal(result, expected[['Y', 'X']])

    @pytest.mark.parametrize('df1_vals, df2_vals', [
        ([0, 1, 2], ["0", "1", "2"]),
        ([0.0, 1.0, 2.0], ["0", "1", "2"]),
        ([0, 1, 2], [u"0", u"1", u"2"]),
        (pd.date_range('1/1/2011', periods=2, freq='D'), ['2011-01-01',
                                                          '2011-01-02']),
        (pd.date_range('1/1/2011', periods=2, freq='D'), [0, 1]),
        (pd.date_range('1/1/2011', periods=2, freq='D'), [0.0, 1.0]),
        (pd.date_range('20130101', periods=3),
            pd.date_range('20130101', periods=3, tz='US/Eastern')),
        ([0, 1, 2], Series(['a', 'b', 'a']).astype('category')),
        ([0.0, 1.0, 2.0], Series(['a', 'b', 'a']).astype('category')),
    ])
    def test_merge_incompat_dtypes(self, df1_vals, df2_vals):
        # GH 9780, GH 15800
        # Raise a ValueError when a user tries to merge on
        # dtypes that are incompatible (e.g., obj and int/float)

        df1 = DataFrame({'A': df1_vals})
        df2 = DataFrame({'A': df2_vals})

        msg = ("You are trying to merge on {lk_dtype} and "
               "{rk_dtype} columns. If you wish to proceed "
               "you should use pd.concat".format(lk_dtype=df1['A'].dtype,
                                                 rk_dtype=df2['A'].dtype))
        msg = re.escape(msg)
        with tm.assert_raises_regex(ValueError, msg):
            pd.merge(df1, df2, on=['A'])

        # Check that error still raised when swapping order of dataframes
        msg = ("You are trying to merge on {lk_dtype} and "
               "{rk_dtype} columns. If you wish to proceed "
               "you should use pd.concat".format(lk_dtype=df2['A'].dtype,
                                                 rk_dtype=df1['A'].dtype))
        msg = re.escape(msg)
        with tm.assert_raises_regex(ValueError, msg):
            pd.merge(df2, df1, on=['A'])


@pytest.fixture
def left():
    np.random.seed(1234)
    return DataFrame(
        {'X': Series(np.random.choice(
            ['foo', 'bar'],
            size=(10,))).astype(CDT(['foo', 'bar'])),
         'Y': np.random.choice(['one', 'two', 'three'], size=(10,))})


@pytest.fixture
def right():
    np.random.seed(1234)
    return DataFrame(
        {'X': Series(['foo', 'bar']).astype(CDT(['foo', 'bar'])),
         'Z': [1, 2]})


class TestMergeCategorical(object):

    def test_identical(self, left):
        # merging on the same, should preserve dtypes
        merged = pd.merge(left, left, on='X')
        result = merged.dtypes.sort_index()
        expected = Series([CategoricalDtype(),
                           np.dtype('O'),
                           np.dtype('O')],
                          index=['X', 'Y_x', 'Y_y'])
        assert_series_equal(result, expected)

    def test_basic(self, left, right):
        # we have matching Categorical dtypes in X
        # so should preserve the merged column
        merged = pd.merge(left, right, on='X')
        result = merged.dtypes.sort_index()
        expected = Series([CategoricalDtype(),
                           np.dtype('O'),
                           np.dtype('int64')],
                          index=['X', 'Y', 'Z'])
        assert_series_equal(result, expected)

    def test_merge_categorical(self):
        # GH 9426

        right = DataFrame({'c': {0: 'a',
                                 1: 'b',
                                 2: 'c',
                                 3: 'd',
                                 4: 'e'},
                           'd': {0: 'null',
                                 1: 'null',
                                 2: 'null',
                                 3: 'null',
                                 4: 'null'}})
        left = DataFrame({'a': {0: 'f',
                                1: 'f',
                                2: 'f',
                                3: 'f',
                                4: 'f'},
                          'b': {0: 'g',
                                1: 'g',
                                2: 'g',
                                3: 'g',
                                4: 'g'}})
        df = pd.merge(left, right, how='left', left_on='b', right_on='c')

        # object-object
        expected = df.copy()

        # object-cat
        # note that we propagate the category
        # because we don't have any matching rows
        cright = right.copy()
        cright['d'] = cright['d'].astype('category')
        result = pd.merge(left, cright, how='left', left_on='b', right_on='c')
        expected['d'] = expected['d'].astype(CategoricalDtype(['null']))
        tm.assert_frame_equal(result, expected)

        # cat-object
        cleft = left.copy()
        cleft['b'] = cleft['b'].astype('category')
        result = pd.merge(cleft, cright, how='left', left_on='b', right_on='c')
        tm.assert_frame_equal(result, expected)

        # cat-cat
        cright = right.copy()
        cright['d'] = cright['d'].astype('category')
        cleft = left.copy()
        cleft['b'] = cleft['b'].astype('category')
        result = pd.merge(cleft, cright, how='left', left_on='b', right_on='c')
        tm.assert_frame_equal(result, expected)

    def tests_merge_categorical_unordered_equal(self):
        # GH-19551
        df1 = DataFrame({
            'Foo': Categorical(['A', 'B', 'C'], categories=['A', 'B', 'C']),
            'Left': ['A0', 'B0', 'C0'],
        })

        df2 = DataFrame({
            'Foo': Categorical(['C', 'B', 'A'], categories=['C', 'B', 'A']),
            'Right': ['C1', 'B1', 'A1'],
        })
        result = pd.merge(df1, df2, on=['Foo'])
        expected = DataFrame({
            'Foo': pd.Categorical(['A', 'B', 'C']),
            'Left': ['A0', 'B0', 'C0'],
            'Right': ['A1', 'B1', 'C1'],
        })
        assert_frame_equal(result, expected)

    def test_other_columns(self, left, right):
        # non-merge columns should preserve if possible
        right = right.assign(Z=right.Z.astype('category'))

        merged = pd.merge(left, right, on='X')
        result = merged.dtypes.sort_index()
        expected = Series([CategoricalDtype(),
                           np.dtype('O'),
                           CategoricalDtype()],
                          index=['X', 'Y', 'Z'])
        assert_series_equal(result, expected)

        # categories are preserved
        assert left.X.values.is_dtype_equal(merged.X.values)
        assert right.Z.values.is_dtype_equal(merged.Z.values)

    @pytest.mark.parametrize(
        'change', [lambda x: x,
                   lambda x: x.astype(CDT(['foo', 'bar', 'bah'])),
                   lambda x: x.astype(CDT(ordered=True))])
    def test_dtype_on_merged_different(self, change, join_type, left, right):
        # our merging columns, X now has 2 different dtypes
        # so we must be object as a result

        X = change(right.X.astype('object'))
        right = right.assign(X=X)
        assert is_categorical_dtype(left.X.values)
        # assert not left.X.values.is_dtype_equal(right.X.values)

        merged = pd.merge(left, right, on='X', how=join_type)

        result = merged.dtypes.sort_index()
        expected = Series([np.dtype('O'),
                           np.dtype('O'),
                           np.dtype('int64')],
                          index=['X', 'Y', 'Z'])
        assert_series_equal(result, expected)

    def test_self_join_multiple_categories(self):
        # GH 16767
        # non-duplicates should work with multiple categories
        m = 5
        df = pd.DataFrame({
            'a': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'] * m,
            'b': ['t', 'w', 'x', 'y', 'z'] * 2 * m,
            'c': [letter
                  for each in ['m', 'n', 'u', 'p', 'o']
                  for letter in [each] * 2 * m],
            'd': [letter
                  for each in ['aa', 'bb', 'cc', 'dd', 'ee',
                               'ff', 'gg', 'hh', 'ii', 'jj']
                  for letter in [each] * m]})

        # change them all to categorical variables
        df = df.apply(lambda x: x.astype('category'))

        # self-join should equal ourselves
        result = pd.merge(df, df, on=list(df.columns))

        assert_frame_equal(result, df)

    def test_dtype_on_categorical_dates(self):
        # GH 16900
        # dates should not be coerced to ints

        df = pd.DataFrame(
            [[date(2001, 1, 1), 1.1],
             [date(2001, 1, 2), 1.3]],
            columns=['date', 'num2']
        )
        df['date'] = df['date'].astype('category')

        df2 = pd.DataFrame(
            [[date(2001, 1, 1), 1.3],
             [date(2001, 1, 3), 1.4]],
            columns=['date', 'num4']
        )
        df2['date'] = df2['date'].astype('category')

        expected_outer = pd.DataFrame([
            [pd.Timestamp('2001-01-01'), 1.1, 1.3],
            [pd.Timestamp('2001-01-02'), 1.3, np.nan],
            [pd.Timestamp('2001-01-03'), np.nan, 1.4]],
            columns=['date', 'num2', 'num4']
        )
        result_outer = pd.merge(df, df2, how='outer', on=['date'])
        assert_frame_equal(result_outer, expected_outer)

        expected_inner = pd.DataFrame(
            [[pd.Timestamp('2001-01-01'), 1.1, 1.3]],
            columns=['date', 'num2', 'num4']
        )
        result_inner = pd.merge(df, df2, how='inner', on=['date'])
        assert_frame_equal(result_inner, expected_inner)

    @pytest.mark.parametrize('ordered', [True, False])
    @pytest.mark.parametrize('category_column,categories,expected_categories',
                             [([False, True, True, False], [True, False],
                               [True, False]),
                              ([2, 1, 1, 2], [1, 2], [1, 2]),
                              (['False', 'True', 'True', 'False'],
                               ['True', 'False'], ['True', 'False'])])
    def test_merging_with_bool_or_int_cateorical_column(self, category_column,
                                                        categories,
                                                        expected_categories,
                                                        ordered):
        # GH 17187
        # merging with a boolean/int categorical column
        df1 = pd.DataFrame({'id': [1, 2, 3, 4],
                            'cat': category_column})
        df1['cat'] = df1['cat'].astype(CDT(categories, ordered=ordered))
        df2 = pd.DataFrame({'id': [2, 4], 'num': [1, 9]})
        result = df1.merge(df2)
        expected = pd.DataFrame({'id': [2, 4], 'cat': expected_categories,
                                 'num': [1, 9]})
        expected['cat'] = expected['cat'].astype(
            CDT(categories, ordered=ordered))
        assert_frame_equal(expected, result)


@pytest.fixture
def left_df():
    return DataFrame({'a': [20, 10, 0]}, index=[2, 1, 0])


@pytest.fixture
def right_df():
    return DataFrame({'b': [300, 100, 200]}, index=[3, 1, 2])


class TestMergeOnIndexes(object):

    @pytest.mark.parametrize(
        "how, sort, expected",
        [('inner', False, DataFrame({'a': [20, 10],
                                     'b': [200, 100]},
                                    index=[2, 1])),
         ('inner', True, DataFrame({'a': [10, 20],
                                    'b': [100, 200]},
                                   index=[1, 2])),
         ('left', False, DataFrame({'a': [20, 10, 0],
                                    'b': [200, 100, np.nan]},
                                   index=[2, 1, 0])),
         ('left', True, DataFrame({'a': [0, 10, 20],
                                   'b': [np.nan, 100, 200]},
                                  index=[0, 1, 2])),
         ('right', False, DataFrame({'a': [np.nan, 10, 20],
                                     'b': [300, 100, 200]},
                                    index=[3, 1, 2])),
         ('right', True, DataFrame({'a': [10, 20, np.nan],
                                    'b': [100, 200, 300]},
                                   index=[1, 2, 3])),
         ('outer', False, DataFrame({'a': [0, 10, 20, np.nan],
                                     'b': [np.nan, 100, 200, 300]},
                                    index=[0, 1, 2, 3])),
         ('outer', True, DataFrame({'a': [0, 10, 20, np.nan],
                                    'b': [np.nan, 100, 200, 300]},
                                   index=[0, 1, 2, 3]))])
    def test_merge_on_indexes(self, left_df, right_df, how, sort, expected):
        result = pd.merge(left_df, right_df,
                          left_index=True,
                          right_index=True,
                          how=how,
                          sort=sort)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    'index', [
        CategoricalIndex(['A', 'B'], categories=['A', 'B'], name='index_col'),
        Float64Index([1.0, 2.0], name='index_col'),
        Int64Index([1, 2], name='index_col'),
        UInt64Index([1, 2], name='index_col'),
        RangeIndex(start=0, stop=2, name='index_col'),
        DatetimeIndex(["2018-01-01", "2018-01-02"], name='index_col'),
    ], ids=lambda x: type(x).__name__)
def test_merge_index_types(index):
    # gh-20777
    # assert key access is consistent across index types
    left = DataFrame({"left_data": [1, 2]}, index=index)
    right = DataFrame({"right_data": [1.0, 2.0]}, index=index)

    result = left.merge(right, on=['index_col'])

    expected = DataFrame(
        OrderedDict([('left_data', [1, 2]), ('right_data', [1.0, 2.0])]),
        index=index)
    assert_frame_equal(result, expected)
