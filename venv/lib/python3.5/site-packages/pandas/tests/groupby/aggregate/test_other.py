# -*- coding: utf-8 -*-

"""
test all other .agg behavior
"""

from __future__ import print_function

import pytest
from collections import OrderedDict

import datetime as dt
from functools import partial

import numpy as np
import pandas as pd

from pandas import (
    date_range, DataFrame, Index, MultiIndex, PeriodIndex, period_range, Series
)
from pandas.core.groupby.groupby import SpecificationError
from pandas.io.formats.printing import pprint_thing
import pandas.util.testing as tm


def test_agg_api():
    # GH 6337
    # http://stackoverflow.com/questions/21706030/pandas-groupby-agg-function-column-dtype-error
    # different api for agg when passed custom function with mixed frame

    df = DataFrame({'data1': np.random.randn(5),
                    'data2': np.random.randn(5),
                    'key1': ['a', 'a', 'b', 'b', 'a'],
                    'key2': ['one', 'two', 'one', 'two', 'one']})
    grouped = df.groupby('key1')

    def peak_to_peak(arr):
        return arr.max() - arr.min()

    expected = grouped.agg([peak_to_peak])
    expected.columns = ['data1', 'data2']
    result = grouped.agg(peak_to_peak)
    tm.assert_frame_equal(result, expected)


def test_agg_datetimes_mixed():
    data = [[1, '2012-01-01', 1.0],
            [2, '2012-01-02', 2.0],
            [3, None, 3.0]]

    df1 = DataFrame({'key': [x[0] for x in data],
                     'date': [x[1] for x in data],
                     'value': [x[2] for x in data]})

    data = [[row[0],
             (dt.datetime.strptime(row[1], '%Y-%m-%d').date()
              if row[1] else None),
             row[2]]
            for row in data]

    df2 = DataFrame({'key': [x[0] for x in data],
                     'date': [x[1] for x in data],
                     'value': [x[2] for x in data]})

    df1['weights'] = df1['value'] / df1['value'].sum()
    gb1 = df1.groupby('date').aggregate(np.sum)

    df2['weights'] = df1['value'] / df1['value'].sum()
    gb2 = df2.groupby('date').aggregate(np.sum)

    assert (len(gb1) == len(gb2))


def test_agg_period_index():
    prng = period_range('2012-1-1', freq='M', periods=3)
    df = DataFrame(np.random.randn(3, 2), index=prng)
    rs = df.groupby(level=0).sum()
    assert isinstance(rs.index, PeriodIndex)

    # GH 3579
    index = period_range(start='1999-01', periods=5, freq='M')
    s1 = Series(np.random.rand(len(index)), index=index)
    s2 = Series(np.random.rand(len(index)), index=index)
    series = [('s1', s1), ('s2', s2)]
    df = DataFrame.from_dict(OrderedDict(series))
    grouped = df.groupby(df.index.month)
    list(grouped)


def test_agg_dict_parameter_cast_result_dtypes():
    # GH 12821

    df = DataFrame({'class': ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'],
                    'time': date_range('1/1/2011', periods=8, freq='H')})
    df.loc[[0, 1, 2, 5], 'time'] = None

    # test for `first` function
    exp = df.loc[[0, 3, 4, 6]].set_index('class')
    grouped = df.groupby('class')
    tm.assert_frame_equal(grouped.first(), exp)
    tm.assert_frame_equal(grouped.agg('first'), exp)
    tm.assert_frame_equal(grouped.agg({'time': 'first'}), exp)
    tm.assert_series_equal(grouped.time.first(), exp['time'])
    tm.assert_series_equal(grouped.time.agg('first'), exp['time'])

    # test for `last` function
    exp = df.loc[[0, 3, 4, 7]].set_index('class')
    grouped = df.groupby('class')
    tm.assert_frame_equal(grouped.last(), exp)
    tm.assert_frame_equal(grouped.agg('last'), exp)
    tm.assert_frame_equal(grouped.agg({'time': 'last'}), exp)
    tm.assert_series_equal(grouped.time.last(), exp['time'])
    tm.assert_series_equal(grouped.time.agg('last'), exp['time'])

    # count
    exp = pd.Series([2, 2, 2, 2],
                    index=Index(list('ABCD'), name='class'),
                    name='time')
    tm.assert_series_equal(grouped.time.agg(len), exp)
    tm.assert_series_equal(grouped.time.size(), exp)

    exp = pd.Series([0, 1, 1, 2],
                    index=Index(list('ABCD'), name='class'),
                    name='time')
    tm.assert_series_equal(grouped.time.count(), exp)


def test_agg_cast_results_dtypes():
    # similar to GH12821
    # xref #11444
    u = [dt.datetime(2015, x + 1, 1) for x in range(12)]
    v = list('aaabbbbbbccd')
    df = pd.DataFrame({'X': v, 'Y': u})

    result = df.groupby('X')['Y'].agg(len)
    expected = df.groupby('X')['Y'].count()
    tm.assert_series_equal(result, expected)


def test_aggregate_float64_no_int64():
    # see gh-11199
    df = DataFrame({"a": [1, 2, 3, 4, 5],
                    "b": [1, 2, 2, 4, 5],
                    "c": [1, 2, 3, 4, 5]})

    expected = DataFrame({"a": [1, 2.5, 4, 5]}, index=[1, 2, 4, 5])
    expected.index.name = "b"

    result = df.groupby("b")[["a"]].mean()
    tm.assert_frame_equal(result, expected)

    expected = DataFrame({"a": [1, 2.5, 4, 5], "c": [1, 2.5, 4, 5]},
                         index=[1, 2, 4, 5])
    expected.index.name = "b"

    result = df.groupby("b")[["a", "c"]].mean()
    tm.assert_frame_equal(result, expected)


def test_aggregate_api_consistency():
    # GH 9052
    # make sure that the aggregates via dict
    # are consistent
    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': np.random.randn(8) + 1.0,
                    'D': np.arange(8)})

    grouped = df.groupby(['A', 'B'])
    c_mean = grouped['C'].mean()
    c_sum = grouped['C'].sum()
    d_mean = grouped['D'].mean()
    d_sum = grouped['D'].sum()

    result = grouped['D'].agg(['sum', 'mean'])
    expected = pd.concat([d_sum, d_mean], axis=1)
    expected.columns = ['sum', 'mean']
    tm.assert_frame_equal(result, expected, check_like=True)

    result = grouped.agg([np.sum, np.mean])
    expected = pd.concat([c_sum, c_mean, d_sum, d_mean], axis=1)
    expected.columns = MultiIndex.from_product([['C', 'D'],
                                                ['sum', 'mean']])
    tm.assert_frame_equal(result, expected, check_like=True)

    result = grouped[['D', 'C']].agg([np.sum, np.mean])
    expected = pd.concat([d_sum, d_mean, c_sum, c_mean], axis=1)
    expected.columns = MultiIndex.from_product([['D', 'C'],
                                                ['sum', 'mean']])
    tm.assert_frame_equal(result, expected, check_like=True)

    result = grouped.agg({'C': 'mean', 'D': 'sum'})
    expected = pd.concat([d_sum, c_mean], axis=1)
    tm.assert_frame_equal(result, expected, check_like=True)

    result = grouped.agg({'C': ['mean', 'sum'],
                          'D': ['mean', 'sum']})
    expected = pd.concat([c_mean, c_sum, d_mean, d_sum], axis=1)
    expected.columns = MultiIndex.from_product([['C', 'D'],
                                                ['mean', 'sum']])

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = grouped[['D', 'C']].agg({'r': np.sum,
                                          'r2': np.mean})
    expected = pd.concat([d_sum, c_sum, d_mean, c_mean], axis=1)
    expected.columns = MultiIndex.from_product([['r', 'r2'],
                                                ['D', 'C']])
    tm.assert_frame_equal(result, expected, check_like=True)


def test_agg_dict_renaming_deprecation():
    # 15931
    df = pd.DataFrame({'A': [1, 1, 1, 2, 2],
                       'B': range(5),
                       'C': range(5)})

    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False) as w:
        df.groupby('A').agg({'B': {'foo': ['sum', 'max']},
                             'C': {'bar': ['count', 'min']}})
        assert "using a dict with renaming" in str(w[0].message)

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        df.groupby('A')[['B', 'C']].agg({'ma': 'max'})

    with tm.assert_produces_warning(FutureWarning) as w:
        df.groupby('A').B.agg({'foo': 'count'})
        assert "using a dict on a Series for aggregation" in str(w[0].message)


def test_agg_compat():
    # GH 12334
    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': np.random.randn(8) + 1.0,
                    'D': np.arange(8)})

    g = df.groupby(['A', 'B'])

    expected = pd.concat([g['D'].sum(), g['D'].std()], axis=1)
    expected.columns = MultiIndex.from_tuples([('C', 'sum'),
                                               ('C', 'std')])
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = g['D'].agg({'C': ['sum', 'std']})
    tm.assert_frame_equal(result, expected, check_like=True)

    expected = pd.concat([g['D'].sum(), g['D'].std()], axis=1)
    expected.columns = ['C', 'D']

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = g['D'].agg({'C': 'sum', 'D': 'std'})
    tm.assert_frame_equal(result, expected, check_like=True)


def test_agg_nested_dicts():
    # API change for disallowing these types of nested dicts
    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                    'B': ['one', 'one', 'two', 'two',
                          'two', 'two', 'one', 'two'],
                    'C': np.random.randn(8) + 1.0,
                    'D': np.arange(8)})

    g = df.groupby(['A', 'B'])

    msg = r'cannot perform renaming for r[1-2] with a nested dictionary'
    with tm.assert_raises_regex(SpecificationError, msg):
        g.aggregate({'r1': {'C': ['mean', 'sum']},
                     'r2': {'D': ['mean', 'sum']}})

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = g.agg({'C': {'ra': ['mean', 'std']},
                        'D': {'rb': ['mean', 'std']}})
    expected = pd.concat([g['C'].mean(), g['C'].std(),
                          g['D'].mean(), g['D'].std()],
                         axis=1)
    expected.columns = pd.MultiIndex.from_tuples(
        [('ra', 'mean'), ('ra', 'std'),
         ('rb', 'mean'), ('rb', 'std')])
    tm.assert_frame_equal(result, expected, check_like=True)

    # same name as the original column
    # GH9052
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        expected = g['D'].agg({'result1': np.sum, 'result2': np.mean})
    expected = expected.rename(columns={'result1': 'D'})

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = g['D'].agg({'D': np.sum, 'result2': np.mean})
    tm.assert_frame_equal(result, expected, check_like=True)


def test_agg_item_by_item_raise_typeerror():
    df = DataFrame(np.random.randint(10, size=(20, 10)))

    def raiseException(df):
        pprint_thing('----------------------------------------')
        pprint_thing(df.to_string())
        raise TypeError('test')

    with tm.assert_raises_regex(TypeError, 'test'):
        df.groupby(0).agg(raiseException)


def test_series_agg_multikey():
    ts = tm.makeTimeSeries()
    grouped = ts.groupby([lambda x: x.year, lambda x: x.month])

    result = grouped.agg(np.sum)
    expected = grouped.sum()
    tm.assert_series_equal(result, expected)


def test_series_agg_multi_pure_python():
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

    def bad(x):
        assert (len(x.values.base) > 0)
        return 'foo'

    result = data.groupby(['A', 'B']).agg(bad)
    expected = data.groupby(['A', 'B']).agg(lambda x: 'foo')
    tm.assert_frame_equal(result, expected)


def test_agg_consistency():
    # agg with ([]) and () not consistent
    # GH 6715
    def P1(a):
        try:
            return np.percentile(a.dropna(), q=1)
        except Exception:
            return np.nan

    df = DataFrame({'col1': [1, 2, 3, 4],
                    'col2': [10, 25, 26, 31],
                    'date': [dt.date(2013, 2, 10), dt.date(2013, 2, 10),
                             dt.date(2013, 2, 11), dt.date(2013, 2, 11)]})

    g = df.groupby('date')

    expected = g.agg([P1])
    expected.columns = expected.columns.levels[0]

    result = g.agg(P1)
    tm.assert_frame_equal(result, expected)


def test_agg_callables():
    # GH 7929
    df = DataFrame({'foo': [1, 2], 'bar': [3, 4]}).astype(np.int64)

    class fn_class(object):

        def __call__(self, x):
            return sum(x)

    equiv_callables = [sum,
                       np.sum,
                       lambda x: sum(x),
                       lambda x: x.sum(),
                       partial(sum),
                       fn_class(), ]

    expected = df.groupby("foo").agg(sum)
    for ecall in equiv_callables:
        result = df.groupby('foo').agg(ecall)
        tm.assert_frame_equal(result, expected)


def test_agg_over_numpy_arrays():
    # GH 3788
    df = pd.DataFrame([[1, np.array([10, 20, 30])],
                       [1, np.array([40, 50, 60])],
                       [2, np.array([20, 30, 40])]],
                      columns=['category', 'arraydata'])
    result = df.groupby('category').agg(sum)

    expected_data = [[np.array([50, 70, 90])], [np.array([20, 30, 40])]]
    expected_index = pd.Index([1, 2], name='category')
    expected_column = ['arraydata']
    expected = pd.DataFrame(expected_data,
                            index=expected_index,
                            columns=expected_column)

    tm.assert_frame_equal(result, expected)


def test_agg_timezone_round_trip():
    # GH 15426
    ts = pd.Timestamp("2016-01-01 12:00:00", tz='US/Pacific')
    df = pd.DataFrame({'a': 1,
                       'b': [ts + dt.timedelta(minutes=nn)
                             for nn in range(10)]})

    result1 = df.groupby('a')['b'].agg(np.min).iloc[0]
    result2 = df.groupby('a')['b'].agg(lambda x: np.min(x)).iloc[0]
    result3 = df.groupby('a')['b'].min().iloc[0]

    assert result1 == ts
    assert result2 == ts
    assert result3 == ts

    dates = [pd.Timestamp("2016-01-0%d 12:00:00" % i, tz='US/Pacific')
             for i in range(1, 5)]
    df = pd.DataFrame({'A': ['a', 'b'] * 2, 'B': dates})
    grouped = df.groupby('A')

    ts = df['B'].iloc[0]
    assert ts == grouped.nth(0)['B'].iloc[0]
    assert ts == grouped.head(1)['B'].iloc[0]
    assert ts == grouped.first()['B'].iloc[0]
    assert ts == grouped.apply(lambda x: x.iloc[0])[0]

    ts = df['B'].iloc[2]
    assert ts == grouped.last()['B'].iloc[0]
    assert ts == grouped.apply(lambda x: x.iloc[-1])[0]


def test_sum_uint64_overflow():
    # see gh-14758
    # Convert to uint64 and don't overflow
    df = pd.DataFrame([[1, 2], [3, 4], [5, 6]], dtype=object)
    df = df + 9223372036854775807

    index = pd.Index([9223372036854775808,
                      9223372036854775810,
                      9223372036854775812],
                     dtype=np.uint64)
    expected = pd.DataFrame({1: [9223372036854775809,
                                 9223372036854775811,
                                 9223372036854775813]},
                            index=index)

    expected.index.name = 0
    result = df.groupby(0).sum()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("structure, expected", [
    (tuple, pd.DataFrame({'C': {(1, 1): (1, 1, 1), (3, 4): (3, 4, 4)}})),
    (list, pd.DataFrame({'C': {(1, 1): [1, 1, 1], (3, 4): [3, 4, 4]}})),
    (lambda x: tuple(x), pd.DataFrame({'C': {(1, 1): (1, 1, 1),
                                             (3, 4): (3, 4, 4)}})),
    (lambda x: list(x), pd.DataFrame({'C': {(1, 1): [1, 1, 1],
                                            (3, 4): [3, 4, 4]}}))
])
def test_agg_structs_dataframe(structure, expected):
    df = pd.DataFrame({'A': [1, 1, 1, 3, 3, 3],
                       'B': [1, 1, 1, 4, 4, 4],
                       'C': [1, 1, 1, 3, 4, 4]})

    result = df.groupby(['A', 'B']).aggregate(structure)
    expected.index.names = ['A', 'B']
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("structure, expected", [
    (tuple, pd.Series([(1, 1, 1), (3, 4, 4)], index=[1, 3], name='C')),
    (list, pd.Series([[1, 1, 1], [3, 4, 4]], index=[1, 3], name='C')),
    (lambda x: tuple(x), pd.Series([(1, 1, 1), (3, 4, 4)],
                                   index=[1, 3], name='C')),
    (lambda x: list(x), pd.Series([[1, 1, 1], [3, 4, 4]],
                                  index=[1, 3], name='C'))
])
def test_agg_structs_series(structure, expected):
    # Issue #18079
    df = pd.DataFrame({'A': [1, 1, 1, 3, 3, 3],
                       'B': [1, 1, 1, 4, 4, 4],
                       'C': [1, 1, 1, 3, 4, 4]})

    result = df.groupby('A')['C'].aggregate(structure)
    expected.index.name = 'A'
    tm.assert_series_equal(result, expected)


@pytest.mark.xfail(reason="GH-18869: agg func not called on empty groups.")
def test_agg_category_nansum(observed):
    categories = ['a', 'b', 'c']
    df = pd.DataFrame({"A": pd.Categorical(['a', 'a', 'b'],
                                           categories=categories),
                       'B': [1, 2, 3]})
    result = df.groupby("A", observed=observed).B.agg(np.nansum)
    expected = pd.Series([3, 3, 0],
                         index=pd.CategoricalIndex(['a', 'b', 'c'],
                                                   categories=categories,
                                                   name='A'),
                         name='B')
    if observed:
        expected = expected[expected != 0]
    tm.assert_series_equal(result, expected)
