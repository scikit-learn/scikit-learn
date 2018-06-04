# -*- coding: utf-8 -*-

"""
test .agg behavior / note that .apply is tested generally in test_groupby.py
"""

import pytest

import numpy as np
import pandas as pd

from pandas import concat, DataFrame, Index, MultiIndex, Series
from pandas.core.groupby.groupby import Grouping, SpecificationError
from pandas.compat import OrderedDict
import pandas.util.testing as tm


def test_agg_regression1(tsframe):
    grouped = tsframe.groupby([lambda x: x.year, lambda x: x.month])
    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


def test_agg_must_agg(df):
    grouped = df.groupby('A')['C']

    msg = "Must produce aggregated value"
    with tm.assert_raises_regex(Exception, msg):
        grouped.agg(lambda x: x.describe())
    with tm.assert_raises_regex(Exception, msg):
        grouped.agg(lambda x: x.index[:2])


def test_agg_ser_multi_key(df):
    # TODO(wesm): unused
    ser = df.C  # noqa

    f = lambda x: x.sum()
    results = df.C.groupby([df.A, df.B]).aggregate(f)
    expected = df.groupby(['A', 'B']).sum()['C']
    tm.assert_series_equal(results, expected)


def test_groupby_aggregation_mixed_dtype():

    # GH 6212
    expected = DataFrame({
        'v1': [5, 5, 7, np.nan, 3, 3, 4, 1],
        'v2': [55, 55, 77, np.nan, 33, 33, 44, 11]},
        index=MultiIndex.from_tuples([(1, 95), (1, 99), (2, 95), (2, 99),
                                      ('big', 'damp'),
                                      ('blue', 'dry'),
                                      ('red', 'red'), ('red', 'wet')],
                                     names=['by1', 'by2']))

    df = DataFrame({
        'v1': [1, 3, 5, 7, 8, 3, 5, np.nan, 4, 5, 7, 9],
        'v2': [11, 33, 55, 77, 88, 33, 55, np.nan, 44, 55, 77, 99],
        'by1': ["red", "blue", 1, 2, np.nan, "big", 1, 2, "red", 1, np.nan,
                12],
        'by2': ["wet", "dry", 99, 95, np.nan, "damp", 95, 99, "red", 99,
                np.nan, np.nan]
    })

    g = df.groupby(['by1', 'by2'])
    result = g[['v1', 'v2']].mean()
    tm.assert_frame_equal(result, expected)


def test_agg_apply_corner(ts, tsframe):
    # nothing to group, all NA
    grouped = ts.groupby(ts * np.nan)
    assert ts.dtype == np.float64

    # groupby float64 values results in Float64Index
    exp = Series([], dtype=np.float64,
                 index=pd.Index([], dtype=np.float64))
    tm.assert_series_equal(grouped.sum(), exp)
    tm.assert_series_equal(grouped.agg(np.sum), exp)
    tm.assert_series_equal(grouped.apply(np.sum), exp,
                           check_index_type=False)

    # DataFrame
    grouped = tsframe.groupby(tsframe['A'] * np.nan)
    exp_df = DataFrame(columns=tsframe.columns, dtype=float,
                       index=pd.Index([], dtype=np.float64))
    tm.assert_frame_equal(grouped.sum(), exp_df, check_names=False)
    tm.assert_frame_equal(grouped.agg(np.sum), exp_df, check_names=False)
    tm.assert_frame_equal(grouped.apply(np.sum), exp_df.iloc[:, :0],
                          check_names=False)


def test_agg_grouping_is_list_tuple(ts):
    df = tm.makeTimeDataFrame()

    grouped = df.groupby(lambda x: x.year)
    grouper = grouped.grouper.groupings[0].grouper
    grouped.grouper.groupings[0] = Grouping(ts.index, list(grouper))

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)

    grouped.grouper.groupings[0] = Grouping(ts.index, tuple(grouper))

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


def test_agg_python_multiindex(mframe):
    grouped = mframe.groupby(['A', 'B'])

    result = grouped.agg(np.mean)
    expected = grouped.mean()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('groupbyfunc', [
    lambda x: x.weekday(),
    [lambda x: x.month, lambda x: x.weekday()],
])
def test_aggregate_str_func(tsframe, groupbyfunc):
    grouped = tsframe.groupby(groupbyfunc)

    # single series
    result = grouped['A'].agg('std')
    expected = grouped['A'].std()
    tm.assert_series_equal(result, expected)

    # group frame by function name
    result = grouped.aggregate('var')
    expected = grouped.var()
    tm.assert_frame_equal(result, expected)

    # group frame by function dict
    result = grouped.agg(OrderedDict([['A', 'var'],
                                      ['B', 'std'],
                                      ['C', 'mean'],
                                      ['D', 'sem']]))
    expected = DataFrame(OrderedDict([['A', grouped['A'].var()],
                                      ['B', grouped['B'].std()],
                                      ['C', grouped['C'].mean()],
                                      ['D', grouped['D'].sem()]]))
    tm.assert_frame_equal(result, expected)


def test_aggregate_item_by_item(df):
    grouped = df.groupby('A')

    aggfun = lambda ser: ser.size
    result = grouped.agg(aggfun)
    foo = (df.A == 'foo').sum()
    bar = (df.A == 'bar').sum()
    K = len(result.columns)

    # GH5782
    # odd comparisons can result here, so cast to make easy
    exp = pd.Series(np.array([foo] * K), index=list('BCD'),
                    dtype=np.float64, name='foo')
    tm.assert_series_equal(result.xs('foo'), exp)

    exp = pd.Series(np.array([bar] * K), index=list('BCD'),
                    dtype=np.float64, name='bar')
    tm.assert_almost_equal(result.xs('bar'), exp)

    def aggfun(ser):
        return ser.size

    result = DataFrame().groupby(df.A).agg(aggfun)
    assert isinstance(result, DataFrame)
    assert len(result) == 0


def test_wrap_agg_out(three_group):
    grouped = three_group.groupby(['A', 'B'])

    def func(ser):
        if ser.dtype == np.object:
            raise TypeError
        else:
            return ser.sum()

    result = grouped.aggregate(func)
    exp_grouped = three_group.loc[:, three_group.columns != 'C']
    expected = exp_grouped.groupby(['A', 'B']).aggregate(func)
    tm.assert_frame_equal(result, expected)


def test_agg_multiple_functions_maintain_order(df):
    # GH #610
    funcs = [('mean', np.mean), ('max', np.max), ('min', np.min)]
    result = df.groupby('A')['C'].agg(funcs)
    exp_cols = Index(['mean', 'max', 'min'])

    tm.assert_index_equal(result.columns, exp_cols)


def test_multiple_functions_tuples_and_non_tuples(df):
    # #1359
    funcs = [('foo', 'mean'), 'std']
    ex_funcs = [('foo', 'mean'), ('std', 'std')]

    result = df.groupby('A')['C'].agg(funcs)
    expected = df.groupby('A')['C'].agg(ex_funcs)
    tm.assert_frame_equal(result, expected)

    result = df.groupby('A').agg(funcs)
    expected = df.groupby('A').agg(ex_funcs)
    tm.assert_frame_equal(result, expected)


def test_agg_multiple_functions_too_many_lambdas(df):
    grouped = df.groupby('A')
    funcs = ['mean', lambda x: x.mean(), lambda x: x.std()]

    msg = 'Function names must be unique, found multiple named <lambda>'
    with tm.assert_raises_regex(SpecificationError, msg):
        grouped.agg(funcs)


def test_more_flexible_frame_multi_function(df):
    grouped = df.groupby('A')

    exmean = grouped.agg(OrderedDict([['C', np.mean], ['D', np.mean]]))
    exstd = grouped.agg(OrderedDict([['C', np.std], ['D', np.std]]))

    expected = concat([exmean, exstd], keys=['mean', 'std'], axis=1)
    expected = expected.swaplevel(0, 1, axis=1).sort_index(level=0, axis=1)

    d = OrderedDict([['C', [np.mean, np.std]], ['D', [np.mean, np.std]]])
    result = grouped.aggregate(d)

    tm.assert_frame_equal(result, expected)

    # be careful
    result = grouped.aggregate(OrderedDict([['C', np.mean],
                                            ['D', [np.mean, np.std]]]))
    expected = grouped.aggregate(OrderedDict([['C', np.mean],
                                              ['D', [np.mean, np.std]]]))
    tm.assert_frame_equal(result, expected)

    def foo(x):
        return np.mean(x)

    def bar(x):
        return np.std(x, ddof=1)

    # this uses column selection & renaming
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        d = OrderedDict([['C', np.mean],
                         ['D', OrderedDict([['foo', np.mean],
                                            ['bar', np.std]])]])
        result = grouped.aggregate(d)

    d = OrderedDict([['C', [np.mean]], ['D', [foo, bar]]])
    expected = grouped.aggregate(d)

    tm.assert_frame_equal(result, expected)


def test_multi_function_flexible_mix(df):
    # GH #1268
    grouped = df.groupby('A')

    # Expected
    d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                     ['D', {'sum': 'sum'}]])
    # this uses column selection & renaming
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        expected = grouped.aggregate(d)

    # Test 1
    d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                     ['D', 'sum']])
    # this uses column selection & renaming
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = grouped.aggregate(d)
    tm.assert_frame_equal(result, expected)

    # Test 2
    d = OrderedDict([['C', OrderedDict([['foo', 'mean'], ['bar', 'std']])],
                     ['D', ['sum']]])
    # this uses column selection & renaming
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = grouped.aggregate(d)
    tm.assert_frame_equal(result, expected)
