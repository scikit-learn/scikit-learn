# -*- coding: utf-8 -*-
from __future__ import print_function

import pytest

from warnings import catch_warnings
from datetime import datetime
from decimal import Decimal

from pandas import (date_range, Timestamp,
                    Index, MultiIndex, DataFrame, Series,
                    Panel, DatetimeIndex, read_csv)
from pandas.errors import PerformanceWarning
from pandas.util.testing import (assert_frame_equal,
                                 assert_series_equal, assert_almost_equal)
from pandas.compat import (range, lrange, StringIO, lmap, lzip, map, zip,
                           OrderedDict)
from pandas import compat
from collections import defaultdict
import pandas.core.common as com
import numpy as np

import pandas.util.testing as tm
import pandas as pd


def test_repr():
    # GH18203
    result = repr(pd.Grouper(key='A', level='B'))
    expected = "Grouper(key='A', level='B', axis=0, sort=False)"
    assert result == expected


@pytest.mark.parametrize('dtype', ['int64', 'int32', 'float64', 'float32'])
def test_basic(dtype):

    data = Series(np.arange(9) // 3, index=np.arange(9), dtype=dtype)

    index = np.arange(9)
    np.random.shuffle(index)
    data = data.reindex(index)

    grouped = data.groupby(lambda x: x // 3)

    for k, v in grouped:
        assert len(v) == 3

    agged = grouped.aggregate(np.mean)
    assert agged[1] == 1

    assert_series_equal(agged, grouped.agg(np.mean))  # shorthand
    assert_series_equal(agged, grouped.mean())
    assert_series_equal(grouped.agg(np.sum), grouped.sum())

    expected = grouped.apply(lambda x: x * x.sum())
    transformed = grouped.transform(lambda x: x * x.sum())
    assert transformed[7] == 12
    assert_series_equal(transformed, expected)

    value_grouped = data.groupby(data)
    assert_series_equal(value_grouped.aggregate(np.mean), agged,
                        check_index_type=False)

    # complex agg
    agged = grouped.aggregate([np.mean, np.std])

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        agged = grouped.aggregate({'one': np.mean, 'two': np.std})

    group_constants = {0: 10, 1: 20, 2: 30}
    agged = grouped.agg(lambda x: group_constants[x.name] + x.mean())
    assert agged[1] == 21

    # corner cases
    pytest.raises(Exception, grouped.aggregate, lambda x: x * 2)


def test_groupby_nonobject_dtype(mframe, df_mixed_floats):
    key = mframe.index.labels[0]
    grouped = mframe.groupby(key)
    result = grouped.sum()

    expected = mframe.groupby(key.astype('O')).sum()
    assert_frame_equal(result, expected)

    # GH 3911, mixed frame non-conversion
    df = df_mixed_floats.copy()
    df['value'] = lrange(len(df))

    def max_value(group):
        return group.loc[group['value'].idxmax()]

    applied = df.groupby('A').apply(max_value)
    result = applied.get_dtype_counts().sort_values()
    expected = Series({'float64': 2,
                       'int64': 1,
                       'object': 2}).sort_values()
    assert_series_equal(result, expected)


def test_groupby_return_type():

    # GH2893, return a reduced type
    df1 = DataFrame(
        [{"val1": 1, "val2": 20},
         {"val1": 1, "val2": 19},
         {"val1": 2, "val2": 27},
         {"val1": 2, "val2": 12}
         ])

    def func(dataf):
        return dataf["val2"] - dataf["val2"].mean()

    result = df1.groupby("val1", squeeze=True).apply(func)
    assert isinstance(result, Series)

    df2 = DataFrame(
        [{"val1": 1, "val2": 20},
         {"val1": 1, "val2": 19},
         {"val1": 1, "val2": 27},
         {"val1": 1, "val2": 12}
         ])

    def func(dataf):
        return dataf["val2"] - dataf["val2"].mean()

    result = df2.groupby("val1", squeeze=True).apply(func)
    assert isinstance(result, Series)

    # GH3596, return a consistent type (regression in 0.11 from 0.10.1)
    df = DataFrame([[1, 1], [1, 1]], columns=['X', 'Y'])
    result = df.groupby('X', squeeze=False).count()
    assert isinstance(result, DataFrame)

    # GH5592
    # inconcistent return type
    df = DataFrame(dict(A=['Tiger', 'Tiger', 'Tiger', 'Lamb', 'Lamb',
                           'Pony', 'Pony'], B=Series(
                               np.arange(7), dtype='int64'), C=date_range(
                                   '20130101', periods=7)))

    def f(grp):
        return grp.iloc[0]

    expected = df.groupby('A').first()[['B']]
    result = df.groupby('A').apply(f)[['B']]
    assert_frame_equal(result, expected)

    def f(grp):
        if grp.name == 'Tiger':
            return None
        return grp.iloc[0]

    result = df.groupby('A').apply(f)[['B']]
    e = expected.copy()
    e.loc['Tiger'] = np.nan
    assert_frame_equal(result, e)

    def f(grp):
        if grp.name == 'Pony':
            return None
        return grp.iloc[0]

    result = df.groupby('A').apply(f)[['B']]
    e = expected.copy()
    e.loc['Pony'] = np.nan
    assert_frame_equal(result, e)

    # 5592 revisited, with datetimes
    def f(grp):
        if grp.name == 'Pony':
            return None
        return grp.iloc[0]

    result = df.groupby('A').apply(f)[['C']]
    e = df.groupby('A').first()[['C']]
    e.loc['Pony'] = pd.NaT
    assert_frame_equal(result, e)

    # scalar outputs
    def f(grp):
        if grp.name == 'Pony':
            return None
        return grp.iloc[0].loc['C']

    result = df.groupby('A').apply(f)
    e = df.groupby('A').first()['C'].copy()
    e.loc['Pony'] = np.nan
    e.name = None
    assert_series_equal(result, e)


def test_pass_args_kwargs(ts, tsframe):

    def f(x, q=None, axis=0):
        return np.percentile(x, q, axis=axis)

    g = lambda x: np.percentile(x, 80, axis=0)

    # Series
    ts_grouped = ts.groupby(lambda x: x.month)
    agg_result = ts_grouped.agg(np.percentile, 80, axis=0)
    apply_result = ts_grouped.apply(np.percentile, 80, axis=0)
    trans_result = ts_grouped.transform(np.percentile, 80, axis=0)

    agg_expected = ts_grouped.quantile(.8)
    trans_expected = ts_grouped.transform(g)

    assert_series_equal(apply_result, agg_expected)
    assert_series_equal(agg_result, agg_expected, check_names=False)
    assert_series_equal(trans_result, trans_expected)

    agg_result = ts_grouped.agg(f, q=80)
    apply_result = ts_grouped.apply(f, q=80)
    trans_result = ts_grouped.transform(f, q=80)
    assert_series_equal(agg_result, agg_expected)
    assert_series_equal(apply_result, agg_expected)
    assert_series_equal(trans_result, trans_expected)

    # DataFrame
    df_grouped = tsframe.groupby(lambda x: x.month)
    agg_result = df_grouped.agg(np.percentile, 80, axis=0)
    apply_result = df_grouped.apply(DataFrame.quantile, .8)
    expected = df_grouped.quantile(.8)
    assert_frame_equal(apply_result, expected)
    assert_frame_equal(agg_result, expected, check_names=False)

    agg_result = df_grouped.agg(f, q=80)
    apply_result = df_grouped.apply(DataFrame.quantile, q=.8)
    assert_frame_equal(agg_result, expected, check_names=False)
    assert_frame_equal(apply_result, expected)


def test_len():
    df = tm.makeTimeDataFrame()
    grouped = df.groupby([lambda x: x.year, lambda x: x.month,
                          lambda x: x.day])
    assert len(grouped) == len(df)

    grouped = df.groupby([lambda x: x.year, lambda x: x.month])
    expected = len({(x.year, x.month) for x in df.index})
    assert len(grouped) == expected

    # issue 11016
    df = pd.DataFrame(dict(a=[np.nan] * 3, b=[1, 2, 3]))
    assert len(df.groupby(('a'))) == 0
    assert len(df.groupby(('b'))) == 3
    assert len(df.groupby(['a', 'b'])) == 3


def test_basic_regression():
    # regression
    T = [1.0 * x for x in lrange(1, 10) * 10][:1095]
    result = Series(T, lrange(0, len(T)))

    groupings = np.random.random((1100, ))
    groupings = Series(groupings, lrange(0, len(groupings))) * 10.

    grouped = result.groupby(groupings)
    grouped.mean()


@pytest.mark.parametrize('dtype', ['float64', 'float32', 'int64',
                                   'int32', 'int16', 'int8'])
def test_with_na_groups(dtype):
    index = Index(np.arange(10))
    values = Series(np.ones(10), index, dtype=dtype)
    labels = Series([np.nan, 'foo', 'bar', 'bar', np.nan, np.nan,
                     'bar', 'bar', np.nan, 'foo'], index=index)

    # this SHOULD be an int
    grouped = values.groupby(labels)
    agged = grouped.agg(len)
    expected = Series([4, 2], index=['bar', 'foo'])

    assert_series_equal(agged, expected, check_dtype=False)

    # assert issubclass(agged.dtype.type, np.integer)

    # explicitly return a float from my function
    def f(x):
        return float(len(x))

    agged = grouped.agg(f)
    expected = Series([4, 2], index=['bar', 'foo'])

    assert_series_equal(agged, expected, check_dtype=False)
    assert issubclass(agged.dtype.type, np.dtype(dtype).type)


def test_indices_concatenation_order():

    # GH 2808

    def f1(x):
        y = x[(x.b % 2) == 1] ** 2
        if y.empty:
            multiindex = MultiIndex(levels=[[]] * 2, labels=[[]] * 2,
                                    names=['b', 'c'])
            res = DataFrame(None, columns=['a'], index=multiindex)
            return res
        else:
            y = y.set_index(['b', 'c'])
            return y

    def f2(x):
        y = x[(x.b % 2) == 1] ** 2
        if y.empty:
            return DataFrame()
        else:
            y = y.set_index(['b', 'c'])
            return y

    def f3(x):
        y = x[(x.b % 2) == 1] ** 2
        if y.empty:
            multiindex = MultiIndex(levels=[[]] * 2, labels=[[]] * 2,
                                    names=['foo', 'bar'])
            res = DataFrame(None, columns=['a', 'b'], index=multiindex)
            return res
        else:
            return y

    df = DataFrame({'a': [1, 2, 2, 2], 'b': lrange(4), 'c': lrange(5, 9)})

    df2 = DataFrame({'a': [3, 2, 2, 2], 'b': lrange(4), 'c': lrange(5, 9)})

    # correct result
    result1 = df.groupby('a').apply(f1)
    result2 = df2.groupby('a').apply(f1)
    assert_frame_equal(result1, result2)

    # should fail (not the same number of levels)
    pytest.raises(AssertionError, df.groupby('a').apply, f2)
    pytest.raises(AssertionError, df2.groupby('a').apply, f2)

    # should fail (incorrect shape)
    pytest.raises(AssertionError, df.groupby('a').apply, f3)
    pytest.raises(AssertionError, df2.groupby('a').apply, f3)


def test_attr_wrapper(ts):
    grouped = ts.groupby(lambda x: x.weekday())

    result = grouped.std()
    expected = grouped.agg(lambda x: np.std(x, ddof=1))
    assert_series_equal(result, expected)

    # this is pretty cool
    result = grouped.describe()
    expected = {}
    for name, gp in grouped:
        expected[name] = gp.describe()
    expected = DataFrame(expected).T
    assert_frame_equal(result, expected)

    # get attribute
    result = grouped.dtype
    expected = grouped.agg(lambda x: x.dtype)

    # make sure raises error
    pytest.raises(AttributeError, getattr, grouped, 'foo')


def test_frame_groupby(tsframe):
    grouped = tsframe.groupby(lambda x: x.weekday())

    # aggregate
    aggregated = grouped.aggregate(np.mean)
    assert len(aggregated) == 5
    assert len(aggregated.columns) == 4

    # by string
    tscopy = tsframe.copy()
    tscopy['weekday'] = [x.weekday() for x in tscopy.index]
    stragged = tscopy.groupby('weekday').aggregate(np.mean)
    assert_frame_equal(stragged, aggregated, check_names=False)

    # transform
    grouped = tsframe.head(30).groupby(lambda x: x.weekday())
    transformed = grouped.transform(lambda x: x - x.mean())
    assert len(transformed) == 30
    assert len(transformed.columns) == 4

    # transform propagate
    transformed = grouped.transform(lambda x: x.mean())
    for name, group in grouped:
        mean = group.mean()
        for idx in group.index:
            tm.assert_series_equal(transformed.xs(idx), mean,
                                   check_names=False)

    # iterate
    for weekday, group in grouped:
        assert group.index[0].weekday() == weekday

    # groups / group_indices
    groups = grouped.groups
    indices = grouped.indices

    for k, v in compat.iteritems(groups):
        samething = tsframe.index.take(indices[k])
        assert (samething == v).all()


def test_frame_groupby_columns(tsframe):
    mapping = {'A': 0, 'B': 0, 'C': 1, 'D': 1}
    grouped = tsframe.groupby(mapping, axis=1)

    # aggregate
    aggregated = grouped.aggregate(np.mean)
    assert len(aggregated) == len(tsframe)
    assert len(aggregated.columns) == 2

    # transform
    tf = lambda x: x - x.mean()
    groupedT = tsframe.T.groupby(mapping, axis=0)
    assert_frame_equal(groupedT.transform(tf).T, grouped.transform(tf))

    # iterate
    for k, v in grouped:
        assert len(v.columns) == 2


def test_frame_set_name_single(df):
    grouped = df.groupby('A')

    result = grouped.mean()
    assert result.index.name == 'A'

    result = df.groupby('A', as_index=False).mean()
    assert result.index.name != 'A'

    result = grouped.agg(np.mean)
    assert result.index.name == 'A'

    result = grouped.agg({'C': np.mean, 'D': np.std})
    assert result.index.name == 'A'

    result = grouped['C'].mean()
    assert result.index.name == 'A'
    result = grouped['C'].agg(np.mean)
    assert result.index.name == 'A'
    result = grouped['C'].agg([np.mean, np.std])
    assert result.index.name == 'A'

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = grouped['C'].agg({'foo': np.mean, 'bar': np.std})
    assert result.index.name == 'A'


def test_multi_func(df):
    col1 = df['A']
    col2 = df['B']

    grouped = df.groupby([col1.get, col2.get])
    agged = grouped.mean()
    expected = df.groupby(['A', 'B']).mean()

    # TODO groupby get drops names
    assert_frame_equal(agged.loc[:, ['C', 'D']],
                       expected.loc[:, ['C', 'D']],
                       check_names=False)

    # some "groups" with no data
    df = DataFrame({'v1': np.random.randn(6),
                    'v2': np.random.randn(6),
                    'k1': np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                    'k2': np.array(['1', '1', '1', '2', '2', '2'])},
                   index=['one', 'two', 'three', 'four', 'five', 'six'])
    # only verify that it works for now
    grouped = df.groupby(['k1', 'k2'])
    grouped.agg(np.sum)


def test_multi_key_multiple_functions(df):
    grouped = df.groupby(['A', 'B'])['C']

    agged = grouped.agg([np.mean, np.std])
    expected = DataFrame({'mean': grouped.agg(np.mean),
                          'std': grouped.agg(np.std)})
    assert_frame_equal(agged, expected)


def test_frame_multi_key_function_list():
    data = DataFrame(
        {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
               'foo', 'foo', 'foo'],
         'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
               'two', 'two', 'one'],
         'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
               'dull', 'shiny', 'shiny', 'shiny'],
         'D': np.random.randn(11),
         'E': np.random.randn(11),
         'F': np.random.randn(11)})

    grouped = data.groupby(['A', 'B'])
    funcs = [np.mean, np.std]
    agged = grouped.agg(funcs)
    expected = pd.concat([grouped['D'].agg(funcs), grouped['E'].agg(funcs),
                          grouped['F'].agg(funcs)],
                         keys=['D', 'E', 'F'], axis=1)
    assert (isinstance(agged.index, MultiIndex))
    assert (isinstance(expected.index, MultiIndex))
    assert_frame_equal(agged, expected)


@pytest.mark.parametrize('op', [lambda x: x.sum(), lambda x: x.mean()])
def test_groupby_multiple_columns(df, op):
    data = df
    grouped = data.groupby(['A', 'B'])

    with catch_warnings(record=True):
        result1 = op(grouped)

        expected = defaultdict(dict)
        for n1, gp1 in data.groupby('A'):
            for n2, gp2 in gp1.groupby('B'):
                expected[n1][n2] = op(gp2.loc[:, ['C', 'D']])
        expected = dict((k, DataFrame(v))
                        for k, v in compat.iteritems(expected))
        expected = Panel.fromDict(expected).swapaxes(0, 1)
        expected.major_axis.name, expected.minor_axis.name = 'A', 'B'

        # a little bit crude
        for col in ['C', 'D']:
            result_col = op(grouped[col])
            exp = expected[col]
            pivoted = result1[col].unstack()
            pivoted2 = result_col.unstack()
            assert_frame_equal(pivoted.reindex_like(exp), exp)
            assert_frame_equal(pivoted2.reindex_like(exp), exp)

    # test single series works the same
    result = data['C'].groupby([data['A'], data['B']]).mean()
    expected = data.groupby(['A', 'B']).mean()['C']

    assert_series_equal(result, expected)


def test_groupby_as_index_agg(df):
    grouped = df.groupby('A', as_index=False)

    # single-key

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    assert_frame_equal(result, expected)

    result2 = grouped.agg(OrderedDict([['C', np.mean], ['D', np.sum]]))
    expected2 = grouped.mean()
    expected2['D'] = grouped.sum()['D']
    assert_frame_equal(result2, expected2)

    grouped = df.groupby('A', as_index=True)
    expected3 = grouped['C'].sum()
    expected3 = DataFrame(expected3).rename(columns={'C': 'Q'})

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result3 = grouped['C'].agg({'Q': np.sum})
    assert_frame_equal(result3, expected3)

    # multi-key

    grouped = df.groupby(['A', 'B'], as_index=False)

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    assert_frame_equal(result, expected)

    result2 = grouped.agg(OrderedDict([['C', np.mean], ['D', np.sum]]))
    expected2 = grouped.mean()
    expected2['D'] = grouped.sum()['D']
    assert_frame_equal(result2, expected2)

    expected3 = grouped['C'].sum()
    expected3 = DataFrame(expected3).rename(columns={'C': 'Q'})
    result3 = grouped['C'].agg({'Q': np.sum})
    assert_frame_equal(result3, expected3)

    # GH7115 & GH8112 & GH8582
    df = DataFrame(np.random.randint(0, 100, (50, 3)),
                   columns=['jim', 'joe', 'jolie'])
    ts = Series(np.random.randint(5, 10, 50), name='jim')

    gr = df.groupby(ts)
    gr.nth(0)  # invokes set_selection_from_grouper internally
    assert_frame_equal(gr.apply(sum), df.groupby(ts).apply(sum))

    for attr in ['mean', 'max', 'count', 'idxmax', 'cumsum', 'all']:
        gr = df.groupby(ts, as_index=False)
        left = getattr(gr, attr)()

        gr = df.groupby(ts.values, as_index=True)
        right = getattr(gr, attr)().reset_index(drop=True)

        assert_frame_equal(left, right)


def test_as_index_series_return_frame(df):
    grouped = df.groupby('A', as_index=False)
    grouped2 = df.groupby(['A', 'B'], as_index=False)

    result = grouped['C'].agg(np.sum)
    expected = grouped.agg(np.sum).loc[:, ['A', 'C']]
    assert isinstance(result, DataFrame)
    assert_frame_equal(result, expected)

    result2 = grouped2['C'].agg(np.sum)
    expected2 = grouped2.agg(np.sum).loc[:, ['A', 'B', 'C']]
    assert isinstance(result2, DataFrame)
    assert_frame_equal(result2, expected2)

    result = grouped['C'].sum()
    expected = grouped.sum().loc[:, ['A', 'C']]
    assert isinstance(result, DataFrame)
    assert_frame_equal(result, expected)

    result2 = grouped2['C'].sum()
    expected2 = grouped2.sum().loc[:, ['A', 'B', 'C']]
    assert isinstance(result2, DataFrame)
    assert_frame_equal(result2, expected2)

    # corner case
    pytest.raises(Exception, grouped['C'].__getitem__, 'D')


def test_groupby_as_index_cython(df):
    data = df

    # single-key
    grouped = data.groupby('A', as_index=False)
    result = grouped.mean()
    expected = data.groupby(['A']).mean()
    expected.insert(0, 'A', expected.index)
    expected.index = np.arange(len(expected))
    assert_frame_equal(result, expected)

    # multi-key
    grouped = data.groupby(['A', 'B'], as_index=False)
    result = grouped.mean()
    expected = data.groupby(['A', 'B']).mean()

    arrays = lzip(*expected.index.values)
    expected.insert(0, 'A', arrays[0])
    expected.insert(1, 'B', arrays[1])
    expected.index = np.arange(len(expected))
    assert_frame_equal(result, expected)


def test_groupby_as_index_series_scalar(df):
    grouped = df.groupby(['A', 'B'], as_index=False)

    # GH #421

    result = grouped['C'].agg(len)
    expected = grouped.agg(len).loc[:, ['A', 'B', 'C']]
    assert_frame_equal(result, expected)


def test_groupby_as_index_corner(df, ts):
    pytest.raises(TypeError, ts.groupby, lambda x: x.weekday(),
                  as_index=False)

    pytest.raises(ValueError, df.groupby, lambda x: x.lower(),
                  as_index=False, axis=1)


def test_groupby_multiple_key(df):
    df = tm.makeTimeDataFrame()
    grouped = df.groupby([lambda x: x.year, lambda x: x.month,
                          lambda x: x.day])
    agged = grouped.sum()
    assert_almost_equal(df.values, agged.values)

    grouped = df.T.groupby([lambda x: x.year,
                            lambda x: x.month,
                            lambda x: x.day], axis=1)

    agged = grouped.agg(lambda x: x.sum())
    tm.assert_index_equal(agged.index, df.columns)
    assert_almost_equal(df.T.values, agged.values)

    agged = grouped.agg(lambda x: x.sum())
    assert_almost_equal(df.T.values, agged.values)


def test_groupby_multi_corner(df):
    # test that having an all-NA column doesn't mess you up
    df = df.copy()
    df['bad'] = np.nan
    agged = df.groupby(['A', 'B']).mean()

    expected = df.groupby(['A', 'B']).mean()
    expected['bad'] = np.nan

    assert_frame_equal(agged, expected)


def test_omit_nuisance(df):
    grouped = df.groupby('A')

    result = grouped.mean()
    expected = df.loc[:, ['A', 'C', 'D']].groupby('A').mean()
    assert_frame_equal(result, expected)

    agged = grouped.agg(np.mean)
    exp = grouped.mean()
    assert_frame_equal(agged, exp)

    df = df.loc[:, ['A', 'C', 'D']]
    df['E'] = datetime.now()
    grouped = df.groupby('A')
    result = grouped.agg(np.sum)
    expected = grouped.sum()
    assert_frame_equal(result, expected)

    # won't work with axis = 1
    grouped = df.groupby({'A': 0, 'C': 0, 'D': 1, 'E': 1}, axis=1)
    result = pytest.raises(TypeError, grouped.agg,
                           lambda x: x.sum(0, numeric_only=False))


def test_omit_nuisance_python_multiple(three_group):
    grouped = three_group.groupby(['A', 'B'])

    agged = grouped.agg(np.mean)
    exp = grouped.mean()
    assert_frame_equal(agged, exp)


def test_empty_groups_corner(mframe):
    # handle empty groups
    df = DataFrame({'k1': np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                    'k2': np.array(['1', '1', '1', '2', '2', '2']),
                    'k3': ['foo', 'bar'] * 3,
                    'v1': np.random.randn(6),
                    'v2': np.random.randn(6)})

    grouped = df.groupby(['k1', 'k2'])
    result = grouped.agg(np.mean)
    expected = grouped.mean()
    assert_frame_equal(result, expected)

    grouped = mframe[3:5].groupby(level=0)
    agged = grouped.apply(lambda x: x.mean())
    agged_A = grouped['A'].apply(np.mean)
    assert_series_equal(agged['A'], agged_A)
    assert agged.index.name == 'first'


def test_nonsense_func():
    df = DataFrame([0])
    pytest.raises(Exception, df.groupby, lambda x: x + 'foo')


def test_wrap_aggregated_output_multindex(mframe):
    df = mframe.T
    df['baz', 'two'] = 'peekaboo'

    keys = [np.array([0, 0, 1]), np.array([0, 0, 1])]
    agged = df.groupby(keys).agg(np.mean)
    assert isinstance(agged.columns, MultiIndex)

    def aggfun(ser):
        if ser.name == ('foo', 'one'):
            raise TypeError
        else:
            return ser.sum()

    agged2 = df.groupby(keys).aggregate(aggfun)
    assert len(agged2.columns) + 1 == len(df.columns)


def test_groupby_level_apply(mframe):

    result = mframe.groupby(level=0).count()
    assert result.index.name == 'first'
    result = mframe.groupby(level=1).count()
    assert result.index.name == 'second'

    result = mframe['A'].groupby(level=0).count()
    assert result.index.name == 'first'


def test_groupby_level_mapper(mframe):
    deleveled = mframe.reset_index()

    mapper0 = {'foo': 0, 'bar': 0, 'baz': 1, 'qux': 1}
    mapper1 = {'one': 0, 'two': 0, 'three': 1}

    result0 = mframe.groupby(mapper0, level=0).sum()
    result1 = mframe.groupby(mapper1, level=1).sum()

    mapped_level0 = np.array([mapper0.get(x) for x in deleveled['first']])
    mapped_level1 = np.array([mapper1.get(x) for x in deleveled['second']])
    expected0 = mframe.groupby(mapped_level0).sum()
    expected1 = mframe.groupby(mapped_level1).sum()
    expected0.index.name, expected1.index.name = 'first', 'second'

    assert_frame_equal(result0, expected0)
    assert_frame_equal(result1, expected1)


def test_groupby_level_nonmulti():
    # GH 1313, GH 13901
    s = Series([1, 2, 3, 10, 4, 5, 20, 6],
               Index([1, 2, 3, 1, 4, 5, 2, 6], name='foo'))
    expected = Series([11, 22, 3, 4, 5, 6],
                      Index(range(1, 7), name='foo'))

    result = s.groupby(level=0).sum()
    tm.assert_series_equal(result, expected)
    result = s.groupby(level=[0]).sum()
    tm.assert_series_equal(result, expected)
    result = s.groupby(level=-1).sum()
    tm.assert_series_equal(result, expected)
    result = s.groupby(level=[-1]).sum()
    tm.assert_series_equal(result, expected)

    pytest.raises(ValueError, s.groupby, level=1)
    pytest.raises(ValueError, s.groupby, level=-2)
    pytest.raises(ValueError, s.groupby, level=[])
    pytest.raises(ValueError, s.groupby, level=[0, 0])
    pytest.raises(ValueError, s.groupby, level=[0, 1])
    pytest.raises(ValueError, s.groupby, level=[1])


def test_groupby_complex():
    # GH 12902
    a = Series(data=np.arange(4) * (1 + 2j), index=[0, 0, 1, 1])
    expected = Series((1 + 2j, 5 + 10j))

    result = a.groupby(level=0).sum()
    assert_series_equal(result, expected)

    result = a.sum(level=0)
    assert_series_equal(result, expected)


def test_mutate_groups():

    # GH3380

    df = DataFrame({
        'cat1': ['a'] * 8 + ['b'] * 6,
        'cat2': ['c'] * 2 + ['d'] * 2 + ['e'] * 2 + ['f'] * 2 + ['c'] * 2 +
        ['d'] * 2 + ['e'] * 2,
        'cat3': lmap(lambda x: 'g%s' % x, lrange(1, 15)),
        'val': np.random.randint(100, size=14),
    })

    def f_copy(x):
        x = x.copy()
        x['rank'] = x.val.rank(method='min')
        return x.groupby('cat2')['rank'].min()

    def f_no_copy(x):
        x['rank'] = x.val.rank(method='min')
        return x.groupby('cat2')['rank'].min()

    grpby_copy = df.groupby('cat1').apply(f_copy)
    grpby_no_copy = df.groupby('cat1').apply(f_no_copy)
    assert_series_equal(grpby_copy, grpby_no_copy)


def test_no_mutate_but_looks_like():

    # GH 8467
    # first show's mutation indicator
    # second does not, but should yield the same results
    df = DataFrame({'key': [1, 1, 1, 2, 2, 2, 3, 3, 3], 'value': range(9)})

    result1 = df.groupby('key', group_keys=True).apply(lambda x: x[:].key)
    result2 = df.groupby('key', group_keys=True).apply(lambda x: x.key)
    assert_series_equal(result1, result2)


def test_groupby_series_indexed_differently():
    s1 = Series([5.0, -9.0, 4.0, 100., -5., 55., 6.7],
                index=Index(['a', 'b', 'c', 'd', 'e', 'f', 'g']))
    s2 = Series([1.0, 1.0, 4.0, 5.0, 5.0, 7.0],
                index=Index(['a', 'b', 'd', 'f', 'g', 'h']))

    grouped = s1.groupby(s2)
    agged = grouped.mean()
    exp = s1.groupby(s2.reindex(s1.index).get).mean()
    assert_series_equal(agged, exp)


def test_groupby_with_hier_columns():
    tuples = list(zip(*[['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux',
                         'qux'], ['one', 'two', 'one', 'two', 'one', 'two',
                                  'one', 'two']]))
    index = MultiIndex.from_tuples(tuples)
    columns = MultiIndex.from_tuples([('A', 'cat'), ('B', 'dog'), (
        'B', 'cat'), ('A', 'dog')])
    df = DataFrame(np.random.randn(8, 4), index=index, columns=columns)

    result = df.groupby(level=0).mean()
    tm.assert_index_equal(result.columns, columns)

    result = df.groupby(level=0, axis=1).mean()
    tm.assert_index_equal(result.index, df.index)

    result = df.groupby(level=0).agg(np.mean)
    tm.assert_index_equal(result.columns, columns)

    result = df.groupby(level=0).apply(lambda x: x.mean())
    tm.assert_index_equal(result.columns, columns)

    result = df.groupby(level=0, axis=1).agg(lambda x: x.mean(1))
    tm.assert_index_equal(result.columns, Index(['A', 'B']))
    tm.assert_index_equal(result.index, df.index)

    # add a nuisance column
    sorted_columns, _ = columns.sortlevel(0)
    df['A', 'foo'] = 'bar'
    result = df.groupby(level=0).mean()
    tm.assert_index_equal(result.columns, df.columns[:-1])


def test_grouping_ndarray(df):
    grouped = df.groupby(df['A'].values)

    result = grouped.sum()
    expected = df.groupby('A').sum()
    assert_frame_equal(result, expected, check_names=False
                       )  # Note: no names when grouping by value


def test_groupby_wrong_multi_labels():
    data = """index,foo,bar,baz,spam,data
0,foo1,bar1,baz1,spam2,20
1,foo1,bar2,baz1,spam3,30
2,foo2,bar2,baz1,spam2,40
3,foo1,bar1,baz2,spam1,50
4,foo3,bar1,baz2,spam1,60"""

    data = read_csv(StringIO(data), index_col=0)

    grouped = data.groupby(['foo', 'bar', 'baz', 'spam'])

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    assert_frame_equal(result, expected)


def test_groupby_series_with_name(df):
    result = df.groupby(df['A']).mean()
    result2 = df.groupby(df['A'], as_index=False).mean()
    assert result.index.name == 'A'
    assert 'A' in result2

    result = df.groupby([df['A'], df['B']]).mean()
    result2 = df.groupby([df['A'], df['B']],
                         as_index=False).mean()
    assert result.index.names == ('A', 'B')
    assert 'A' in result2
    assert 'B' in result2


def test_seriesgroupby_name_attr(df):
    # GH 6265
    result = df.groupby('A')['C']
    assert result.count().name == 'C'
    assert result.mean().name == 'C'

    testFunc = lambda x: np.sum(x) * 2
    assert result.agg(testFunc).name == 'C'


def test_consistency_name():
    # GH 12363

    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': np.random.randn(8) + 1.0,
                    'D': np.arange(8)})

    expected = df.groupby(['A']).B.count()
    result = df.B.groupby(df.A).count()
    assert_series_equal(result, expected)


def test_groupby_name_propagation(df):
    # GH 6124
    def summarize(df, name=None):
        return Series({'count': 1, 'mean': 2, 'omissions': 3, }, name=name)

    def summarize_random_name(df):
        # Provide a different name for each Series.  In this case, groupby
        # should not attempt to propagate the Series name since they are
        # inconsistent.
        return Series({
            'count': 1,
            'mean': 2,
            'omissions': 3,
        }, name=df.iloc[0]['A'])

    metrics = df.groupby('A').apply(summarize)
    assert metrics.columns.name is None
    metrics = df.groupby('A').apply(summarize, 'metrics')
    assert metrics.columns.name == 'metrics'
    metrics = df.groupby('A').apply(summarize_random_name)
    assert metrics.columns.name is None


def test_groupby_nonstring_columns():
    df = DataFrame([np.arange(10) for x in range(10)])
    grouped = df.groupby(0)
    result = grouped.mean()
    expected = df.groupby(df[0]).mean()
    assert_frame_equal(result, expected)


def test_groupby_mixed_type_columns():
    # GH 13432, unorderable types in py3
    df = DataFrame([[0, 1, 2]], columns=['A', 'B', 0])
    expected = DataFrame([[1, 2]], columns=['B', 0],
                         index=Index([0], name='A'))

    result = df.groupby('A').first()
    tm.assert_frame_equal(result, expected)

    result = df.groupby('A').sum()
    tm.assert_frame_equal(result, expected)


def test_cython_grouper_series_bug_noncontig():
    arr = np.empty((100, 100))
    arr.fill(np.nan)
    obj = Series(arr[:, 0], index=lrange(100))
    inds = np.tile(lrange(10), 10)

    result = obj.groupby(inds).agg(Series.median)
    assert result.isna().all()


def test_series_grouper_noncontig_index():
    index = Index(tm.rands_array(10, 100))

    values = Series(np.random.randn(50), index=index[::2])
    labels = np.random.randint(0, 5, 50)

    # it works!
    grouped = values.groupby(labels)

    # accessing the index elements causes segfault
    f = lambda x: len(set(map(id, x.index)))
    grouped.agg(f)


def test_convert_objects_leave_decimal_alone():

    s = Series(lrange(5))
    labels = np.array(['a', 'b', 'c', 'd', 'e'], dtype='O')

    def convert_fast(x):
        return Decimal(str(x.mean()))

    def convert_force_pure(x):
        # base will be length 0
        assert (len(x.values.base) > 0)
        return Decimal(str(x.mean()))

    grouped = s.groupby(labels)

    result = grouped.agg(convert_fast)
    assert result.dtype == np.object_
    assert isinstance(result[0], Decimal)

    result = grouped.agg(convert_force_pure)
    assert result.dtype == np.object_
    assert isinstance(result[0], Decimal)


def test_groupby_dtype_inference_empty():
    # GH 6733
    df = DataFrame({'x': [], 'range': np.arange(0, dtype='int64')})
    assert df['x'].dtype == np.float64

    result = df.groupby('x').first()
    exp_index = Index([], name='x', dtype=np.float64)
    expected = DataFrame({'range': Series(
        [], index=exp_index, dtype='int64')})
    assert_frame_equal(result, expected, by_blocks=True)


def test_groupby_list_infer_array_like(df):
    result = df.groupby(list(df['A'])).mean()
    expected = df.groupby(df['A']).mean()
    assert_frame_equal(result, expected, check_names=False)

    pytest.raises(Exception, df.groupby, list(df['A'][:-1]))

    # pathological case of ambiguity
    df = DataFrame({'foo': [0, 1],
                    'bar': [3, 4],
                    'val': np.random.randn(2)})

    result = df.groupby(['foo', 'bar']).mean()
    expected = df.groupby([df['foo'], df['bar']]).mean()[['val']]


def test_groupby_keys_same_size_as_index():
    # GH 11185
    freq = 's'
    index = pd.date_range(start=pd.Timestamp('2015-09-29T11:34:44-0700'),
                          periods=2, freq=freq)
    df = pd.DataFrame([['A', 10], ['B', 15]], columns=[
        'metric', 'values'
    ], index=index)
    result = df.groupby([pd.Grouper(level=0, freq=freq), 'metric']).mean()
    expected = df.set_index([df.index, 'metric'])

    assert_frame_equal(result, expected)


def test_groupby_one_row():
    # GH 11741
    df1 = pd.DataFrame(np.random.randn(1, 4), columns=list('ABCD'))
    pytest.raises(KeyError, df1.groupby, 'Z')
    df2 = pd.DataFrame(np.random.randn(2, 4), columns=list('ABCD'))
    pytest.raises(KeyError, df2.groupby, 'Z')


def test_groupby_nat_exclude():
    # GH 6992
    df = pd.DataFrame(
        {'values': np.random.randn(8),
         'dt': [np.nan, pd.Timestamp('2013-01-01'), np.nan, pd.Timestamp(
             '2013-02-01'), np.nan, pd.Timestamp('2013-02-01'), np.nan,
            pd.Timestamp('2013-01-01')],
         'str': [np.nan, 'a', np.nan, 'a', np.nan, 'a', np.nan, 'b']})
    grouped = df.groupby('dt')

    expected = [pd.Index([1, 7]), pd.Index([3, 5])]
    keys = sorted(grouped.groups.keys())
    assert len(keys) == 2
    for k, e in zip(keys, expected):
        # grouped.groups keys are np.datetime64 with system tz
        # not to be affected by tz, only compare values
        tm.assert_index_equal(grouped.groups[k], e)

    # confirm obj is not filtered
    tm.assert_frame_equal(grouped.grouper.groupings[0].obj, df)
    assert grouped.ngroups == 2

    expected = {
        Timestamp('2013-01-01 00:00:00'): np.array([1, 7], dtype=np.int64),
        Timestamp('2013-02-01 00:00:00'): np.array([3, 5], dtype=np.int64)
    }

    for k in grouped.indices:
        tm.assert_numpy_array_equal(grouped.indices[k], expected[k])

    tm.assert_frame_equal(
        grouped.get_group(Timestamp('2013-01-01')), df.iloc[[1, 7]])
    tm.assert_frame_equal(
        grouped.get_group(Timestamp('2013-02-01')), df.iloc[[3, 5]])

    pytest.raises(KeyError, grouped.get_group, pd.NaT)

    nan_df = DataFrame({'nan': [np.nan, np.nan, np.nan],
                        'nat': [pd.NaT, pd.NaT, pd.NaT]})
    assert nan_df['nan'].dtype == 'float64'
    assert nan_df['nat'].dtype == 'datetime64[ns]'

    for key in ['nan', 'nat']:
        grouped = nan_df.groupby(key)
        assert grouped.groups == {}
        assert grouped.ngroups == 0
        assert grouped.indices == {}
        pytest.raises(KeyError, grouped.get_group, np.nan)
        pytest.raises(KeyError, grouped.get_group, pd.NaT)


def test_sparse_friendly(df):
    sdf = df[['C', 'D']].to_sparse()
    with catch_warnings(record=True):
        panel = tm.makePanel()
        tm.add_nans(panel)

    def _check_work(gp):
        gp.mean()
        gp.agg(np.mean)
        dict(iter(gp))

    # it works!
    _check_work(sdf.groupby(lambda x: x // 2))
    _check_work(sdf['C'].groupby(lambda x: x // 2))
    _check_work(sdf.groupby(df['A']))

    # do this someday
    # _check_work(panel.groupby(lambda x: x.month, axis=1))


def test_panel_groupby():
    with catch_warnings(record=True):
        panel = tm.makePanel()
        tm.add_nans(panel)
        grouped = panel.groupby({'ItemA': 0, 'ItemB': 0, 'ItemC': 1},
                                axis='items')
        agged = grouped.mean()
        agged2 = grouped.agg(lambda x: x.mean('items'))

        tm.assert_panel_equal(agged, agged2)

        tm.assert_index_equal(agged.items, Index([0, 1]))

        grouped = panel.groupby(lambda x: x.month, axis='major')
        agged = grouped.mean()

        exp = Index(sorted(list(set(panel.major_axis.month))))
        tm.assert_index_equal(agged.major_axis, exp)

        grouped = panel.groupby({'A': 0, 'B': 0, 'C': 1, 'D': 1},
                                axis='minor')
        agged = grouped.mean()
        tm.assert_index_equal(agged.minor_axis, Index([0, 1]))


def test_groupby_2d_malformed():
    d = DataFrame(index=lrange(2))
    d['group'] = ['g1', 'g2']
    d['zeros'] = [0, 0]
    d['ones'] = [1, 1]
    d['label'] = ['l1', 'l2']
    tmp = d.groupby(['group']).mean()
    res_values = np.array([[0, 1], [0, 1]], dtype=np.int64)
    tm.assert_index_equal(tmp.columns, Index(['zeros', 'ones']))
    tm.assert_numpy_array_equal(tmp.values, res_values)


def test_int32_overflow():
    B = np.concatenate((np.arange(10000), np.arange(10000), np.arange(5000)
                        ))
    A = np.arange(25000)
    df = DataFrame({'A': A,
                    'B': B,
                    'C': A,
                    'D': B,
                    'E': np.random.randn(25000)})

    left = df.groupby(['A', 'B', 'C', 'D']).sum()
    right = df.groupby(['D', 'C', 'B', 'A']).sum()
    assert len(left) == len(right)


def test_groupby_sort_multi():
    df = DataFrame({'a': ['foo', 'bar', 'baz'],
                    'b': [3, 2, 1],
                    'c': [0, 1, 2],
                    'd': np.random.randn(3)})

    tups = lmap(tuple, df[['a', 'b', 'c']].values)
    tups = com._asarray_tuplesafe(tups)
    result = df.groupby(['a', 'b', 'c'], sort=True).sum()
    tm.assert_numpy_array_equal(result.index.values, tups[[1, 2, 0]])

    tups = lmap(tuple, df[['c', 'a', 'b']].values)
    tups = com._asarray_tuplesafe(tups)
    result = df.groupby(['c', 'a', 'b'], sort=True).sum()
    tm.assert_numpy_array_equal(result.index.values, tups)

    tups = lmap(tuple, df[['b', 'c', 'a']].values)
    tups = com._asarray_tuplesafe(tups)
    result = df.groupby(['b', 'c', 'a'], sort=True).sum()
    tm.assert_numpy_array_equal(result.index.values, tups[[2, 1, 0]])

    df = DataFrame({'a': [0, 1, 2, 0, 1, 2],
                    'b': [0, 0, 0, 1, 1, 1],
                    'd': np.random.randn(6)})
    grouped = df.groupby(['a', 'b'])['d']
    result = grouped.sum()

    def _check_groupby(df, result, keys, field, f=lambda x: x.sum()):
        tups = lmap(tuple, df[keys].values)
        tups = com._asarray_tuplesafe(tups)
        expected = f(df.groupby(tups)[field])
        for k, v in compat.iteritems(expected):
            assert (result[k] == v)

    _check_groupby(df, result, ['a', 'b'], 'd')


def test_dont_clobber_name_column():
    df = DataFrame({'key': ['a', 'a', 'a', 'b', 'b', 'b'],
                    'name': ['foo', 'bar', 'baz'] * 2})

    result = df.groupby('key').apply(lambda x: x)
    assert_frame_equal(result, df)


def test_skip_group_keys():

    tsf = tm.makeTimeDataFrame()

    grouped = tsf.groupby(lambda x: x.month, group_keys=False)
    result = grouped.apply(lambda x: x.sort_values(by='A')[:3])

    pieces = []
    for key, group in grouped:
        pieces.append(group.sort_values(by='A')[:3])

    expected = pd.concat(pieces)
    assert_frame_equal(result, expected)

    grouped = tsf['A'].groupby(lambda x: x.month, group_keys=False)
    result = grouped.apply(lambda x: x.sort_values()[:3])

    pieces = []
    for key, group in grouped:
        pieces.append(group.sort_values()[:3])

    expected = pd.concat(pieces)
    assert_series_equal(result, expected)


def test_no_nonsense_name(frame):
    # GH #995
    s = frame['C'].copy()
    s.name = None

    result = s.groupby(frame['A']).agg(np.sum)
    assert result.name is None


def test_multifunc_sum_bug():
    # GH #1065
    x = DataFrame(np.arange(9).reshape(3, 3))
    x['test'] = 0
    x['fl'] = [1.3, 1.5, 1.6]

    grouped = x.groupby('test')
    result = grouped.agg({'fl': 'sum', 2: 'size'})
    assert result['fl'].dtype == np.float64


def test_handle_dict_return_value(df):
    def f(group):
        return {'max': group.max(), 'min': group.min()}

    def g(group):
        return Series({'max': group.max(), 'min': group.min()})

    result = df.groupby('A')['C'].apply(f)
    expected = df.groupby('A')['C'].apply(g)

    assert isinstance(result, Series)
    assert_series_equal(result, expected)


@pytest.mark.parametrize('grouper', ['A', ['A', 'B']])
def test_set_group_name(df, grouper):
    def f(group):
        assert group.name is not None
        return group

    def freduce(group):
        assert group.name is not None
        return group.sum()

    def foo(x):
        return freduce(x)

    grouped = df.groupby(grouper)

    # make sure all these work
    grouped.apply(f)
    grouped.aggregate(freduce)
    grouped.aggregate({'C': freduce, 'D': freduce})
    grouped.transform(f)

    grouped['C'].apply(f)
    grouped['C'].aggregate(freduce)
    grouped['C'].aggregate([freduce, foo])
    grouped['C'].transform(f)


def test_group_name_available_in_inference_pass():
    # gh-15062
    df = pd.DataFrame({'a': [0, 0, 1, 1, 2, 2], 'b': np.arange(6)})

    names = []

    def f(group):
        names.append(group.name)
        return group.copy()

    df.groupby('a', sort=False, group_keys=False).apply(f)
    # we expect 2 zeros because we call ``f`` once to see if a faster route
    # can be used.
    expected_names = [0, 0, 1, 2]
    assert names == expected_names


def test_no_dummy_key_names(df):
    # see gh-1291
    result = df.groupby(df['A'].values).sum()
    assert result.index.name is None

    result = df.groupby([df['A'].values, df['B'].values]).sum()
    assert result.index.names == (None, None)


def test_groupby_sort_multiindex_series():
    # series multiindex groupby sort argument was not being passed through
    # _compress_group_index
    # GH 9444
    index = MultiIndex(levels=[[1, 2], [1, 2]],
                       labels=[[0, 0, 0, 0, 1, 1], [1, 1, 0, 0, 0, 0]],
                       names=['a', 'b'])
    mseries = Series([0, 1, 2, 3, 4, 5], index=index)
    index = MultiIndex(levels=[[1, 2], [1, 2]],
                       labels=[[0, 0, 1], [1, 0, 0]], names=['a', 'b'])
    mseries_result = Series([0, 2, 4], index=index)

    result = mseries.groupby(level=['a', 'b'], sort=False).first()
    assert_series_equal(result, mseries_result)
    result = mseries.groupby(level=['a', 'b'], sort=True).first()
    assert_series_equal(result, mseries_result.sort_index())


def test_groupby_reindex_inside_function():

    periods = 1000
    ind = DatetimeIndex(start='2012/1/1', freq='5min', periods=periods)
    df = DataFrame({'high': np.arange(
        periods), 'low': np.arange(periods)}, index=ind)

    def agg_before(hour, func, fix=False):
        """
            Run an aggregate func on the subset of data.
        """

        def _func(data):
            d = data.loc[data.index.map(
                lambda x: x.hour < 11)].dropna()
            if fix:
                data[data.index[0]]
            if len(d) == 0:
                return None
            return func(d)

        return _func

    def afunc(data):
        d = data.select(lambda x: x.hour < 11).dropna()
        return np.max(d)

    grouped = df.groupby(lambda x: datetime(x.year, x.month, x.day))
    closure_bad = grouped.agg({'high': agg_before(11, np.max)})
    closure_good = grouped.agg({'high': agg_before(11, np.max, True)})

    assert_frame_equal(closure_bad, closure_good)


def test_groupby_multiindex_missing_pair():
    # GH9049
    df = DataFrame({'group1': ['a', 'a', 'a', 'b'],
                    'group2': ['c', 'c', 'd', 'c'],
                    'value': [1, 1, 1, 5]})
    df = df.set_index(['group1', 'group2'])
    df_grouped = df.groupby(level=['group1', 'group2'], sort=True)

    res = df_grouped.agg('sum')
    idx = MultiIndex.from_tuples(
        [('a', 'c'), ('a', 'd'), ('b', 'c')], names=['group1', 'group2'])
    exp = DataFrame([[2], [1], [5]], index=idx, columns=['value'])

    tm.assert_frame_equal(res, exp)


def test_groupby_multiindex_not_lexsorted():
    # GH 11640

    # define the lexsorted version
    lexsorted_mi = MultiIndex.from_tuples(
        [('a', ''), ('b1', 'c1'), ('b2', 'c2')], names=['b', 'c'])
    lexsorted_df = DataFrame([[1, 3, 4]], columns=lexsorted_mi)
    assert lexsorted_df.columns.is_lexsorted()

    # define the non-lexsorted version
    not_lexsorted_df = DataFrame(columns=['a', 'b', 'c', 'd'],
                                 data=[[1, 'b1', 'c1', 3],
                                       [1, 'b2', 'c2', 4]])
    not_lexsorted_df = not_lexsorted_df.pivot_table(
        index='a', columns=['b', 'c'], values='d')
    not_lexsorted_df = not_lexsorted_df.reset_index()
    assert not not_lexsorted_df.columns.is_lexsorted()

    # compare the results
    tm.assert_frame_equal(lexsorted_df, not_lexsorted_df)

    expected = lexsorted_df.groupby('a').mean()
    with tm.assert_produces_warning(PerformanceWarning):
        result = not_lexsorted_df.groupby('a').mean()
    tm.assert_frame_equal(expected, result)

    # a transforming function should work regardless of sort
    # GH 14776
    df = DataFrame({'x': ['a', 'a', 'b', 'a'],
                    'y': [1, 1, 2, 2],
                    'z': [1, 2, 3, 4]}).set_index(['x', 'y'])
    assert not df.index.is_lexsorted()

    for level in [0, 1, [0, 1]]:
        for sort in [False, True]:
            result = df.groupby(level=level, sort=sort).apply(
                DataFrame.drop_duplicates)
            expected = df
            tm.assert_frame_equal(expected, result)

            result = df.sort_index().groupby(level=level, sort=sort).apply(
                DataFrame.drop_duplicates)
            expected = df.sort_index()
            tm.assert_frame_equal(expected, result)


def test_index_label_overlaps_location():
    # checking we don't have any label/location confusion in the
    # the wake of GH5375
    df = DataFrame(list('ABCDE'), index=[2, 0, 2, 1, 1])
    g = df.groupby(list('ababb'))
    actual = g.filter(lambda x: len(x) > 2)
    expected = df.iloc[[1, 3, 4]]
    assert_frame_equal(actual, expected)

    ser = df[0]
    g = ser.groupby(list('ababb'))
    actual = g.filter(lambda x: len(x) > 2)
    expected = ser.take([1, 3, 4])
    assert_series_equal(actual, expected)

    # ... and again, with a generic Index of floats
    df.index = df.index.astype(float)
    g = df.groupby(list('ababb'))
    actual = g.filter(lambda x: len(x) > 2)
    expected = df.iloc[[1, 3, 4]]
    assert_frame_equal(actual, expected)

    ser = df[0]
    g = ser.groupby(list('ababb'))
    actual = g.filter(lambda x: len(x) > 2)
    expected = ser.take([1, 3, 4])
    assert_series_equal(actual, expected)


def test_transform_doesnt_clobber_ints():
    # GH 7972
    n = 6
    x = np.arange(n)
    df = DataFrame({'a': x // 2, 'b': 2.0 * x, 'c': 3.0 * x})
    df2 = DataFrame({'a': x // 2 * 1.0, 'b': 2.0 * x, 'c': 3.0 * x})

    gb = df.groupby('a')
    result = gb.transform('mean')

    gb2 = df2.groupby('a')
    expected = gb2.transform('mean')
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('sort_column', ['ints', 'floats', 'strings',
                                         ['ints', 'floats'],
                                         ['ints', 'strings']])
@pytest.mark.parametrize('group_column', ['int_groups', 'string_groups',
                                          ['int_groups', 'string_groups']])
def test_groupby_preserves_sort(sort_column, group_column):
    # Test to ensure that groupby always preserves sort order of original
    # object. Issue #8588 and #9651

    df = DataFrame(
        {'int_groups': [3, 1, 0, 1, 0, 3, 3, 3],
         'string_groups': ['z', 'a', 'z', 'a', 'a', 'g', 'g', 'g'],
         'ints': [8, 7, 4, 5, 2, 9, 1, 1],
         'floats': [2.3, 5.3, 6.2, -2.4, 2.2, 1.1, 1.1, 5],
         'strings': ['z', 'd', 'a', 'e', 'word', 'word2', '42', '47']})

    # Try sorting on different types and with different group types

    df = df.sort_values(by=sort_column)
    g = df.groupby(group_column)

    def test_sort(x):
        assert_frame_equal(x, x.sort_values(by=sort_column))
    g.apply(test_sort)


def test_group_shift_with_null_key():
    # This test is designed to replicate the segfault in issue #13813.
    n_rows = 1200

    # Generate a moderately large dataframe with occasional missing
    # values in column `B`, and then group by [`A`, `B`]. This should
    # force `-1` in `labels` array of `g.grouper.group_info` exactly
    # at those places, where the group-by key is partially missing.
    df = DataFrame([(i % 12, i % 3 if i % 3 else np.nan, i)
                    for i in range(n_rows)], dtype=float,
                   columns=["A", "B", "Z"], index=None)
    g = df.groupby(["A", "B"])

    expected = DataFrame([(i + 12 if i % 3 and i < n_rows - 12
                           else np.nan)
                          for i in range(n_rows)], dtype=float,
                         columns=["Z"], index=None)
    result = g.shift(-1)

    assert_frame_equal(result, expected)


def test_pivot_table_values_key_error():
    # This test is designed to replicate the error in issue #14938
    df = pd.DataFrame({'eventDate':
                       pd.date_range(pd.datetime.today(),
                                     periods=20, freq='M').tolist(),
                       'thename': range(0, 20)})

    df['year'] = df.set_index('eventDate').index.year
    df['month'] = df.set_index('eventDate').index.month

    with pytest.raises(KeyError):
        df.reset_index().pivot_table(index='year', columns='month',
                                     values='badname', aggfunc='count')


def test_empty_dataframe_groupby():
    # GH8093
    df = DataFrame(columns=['A', 'B', 'C'])

    result = df.groupby('A').sum()
    expected = DataFrame(columns=['B', 'C'], dtype=np.float64)
    expected.index.name = 'A'

    assert_frame_equal(result, expected)


def test_tuple_warns():
    # https://github.com/pandas-dev/pandas/issues/18314
    df = pd.DataFrame({('a', 'b'): [1, 1, 2, 2], 'a': [1, 1, 1, 2],
                       'b': [1, 2, 2, 2], 'c': [1, 1, 1, 1]})
    with tm.assert_produces_warning(FutureWarning) as w:
        df[['a', 'b', 'c']].groupby(('a', 'b')).c.mean()

    assert "Interpreting tuple 'by' as a list" in str(w[0].message)

    with tm.assert_produces_warning(None):
        df.groupby(('a', 'b')).c.mean()


def test_tuple_warns_unhashable():
    # https://github.com/pandas-dev/pandas/issues/18314
    business_dates = date_range(start='4/1/2014', end='6/30/2014',
                                freq='B')
    df = DataFrame(1, index=business_dates, columns=['a', 'b'])

    with tm.assert_produces_warning(FutureWarning) as w:
        df.groupby((df.index.year, df.index.month)).nth([0, 3, -1])

    assert "Interpreting tuple 'by' as a list" in str(w[0].message)


def test_tuple_correct_keyerror():
    # https://github.com/pandas-dev/pandas/issues/18798
    df = pd.DataFrame(1, index=range(3),
                      columns=pd.MultiIndex.from_product([[1, 2],
                                                          [3, 4]]))
    with tm.assert_raises_regex(KeyError, "(7, 8)"):
        df.groupby((7, 8)).mean()
