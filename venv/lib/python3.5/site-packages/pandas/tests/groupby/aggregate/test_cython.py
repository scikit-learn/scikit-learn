# -*- coding: utf-8 -*-

"""
test cython .agg behavior
"""

from __future__ import print_function

import pytest

import numpy as np
from numpy import nan
import pandas as pd

from pandas import (bdate_range, DataFrame, Index, Series, Timestamp,
                    Timedelta, NaT)
from pandas.core.groupby.groupby import DataError
import pandas.util.testing as tm


@pytest.mark.parametrize('op_name', [
    'count',
    'sum',
    'std',
    'var',
    'sem',
    'mean',
    'median',
    'prod',
    'min',
    'max',
])
def test_cythonized_aggers(op_name):
    data = {'A': [0, 0, 0, 0, 1, 1, 1, 1, 1, 1., nan, nan],
            'B': ['A', 'B'] * 6,
            'C': np.random.randn(12)}
    df = DataFrame(data)
    df.loc[2:10:2, 'C'] = nan

    op = lambda x: getattr(x, op_name)()

    # single column
    grouped = df.drop(['B'], axis=1).groupby('A')
    exp = {}
    for cat, group in grouped:
        exp[cat] = op(group['C'])
    exp = DataFrame({'C': exp})
    exp.index.name = 'A'
    result = op(grouped)
    tm.assert_frame_equal(result, exp)

    # multiple columns
    grouped = df.groupby(['A', 'B'])
    expd = {}
    for (cat1, cat2), group in grouped:
        expd.setdefault(cat1, {})[cat2] = op(group['C'])
    exp = DataFrame(expd).T.stack(dropna=False)
    exp.index.names = ['A', 'B']
    exp.name = 'C'

    result = op(grouped)['C']
    if op_name in ['sum', 'prod']:
        tm.assert_series_equal(result, exp)


def test_cython_agg_boolean():
    frame = DataFrame({'a': np.random.randint(0, 5, 50),
                       'b': np.random.randint(0, 2, 50).astype('bool')})
    result = frame.groupby('a')['b'].mean()
    expected = frame.groupby('a')['b'].agg(np.mean)

    tm.assert_series_equal(result, expected)


def test_cython_agg_nothing_to_agg():
    frame = DataFrame({'a': np.random.randint(0, 5, 50),
                       'b': ['foo', 'bar'] * 25})
    msg = "No numeric types to aggregate"

    with tm.assert_raises_regex(DataError, msg):
        frame.groupby('a')['b'].mean()

    frame = DataFrame({'a': np.random.randint(0, 5, 50),
                       'b': ['foo', 'bar'] * 25})
    with tm.assert_raises_regex(DataError, msg):
        frame[['b']].groupby(frame['a']).mean()


def test_cython_agg_nothing_to_agg_with_dates():
    frame = DataFrame({'a': np.random.randint(0, 5, 50),
                       'b': ['foo', 'bar'] * 25,
                       'dates': pd.date_range('now', periods=50, freq='T')})
    msg = "No numeric types to aggregate"
    with tm.assert_raises_regex(DataError, msg):
        frame.groupby('b').dates.mean()


def test_cython_agg_frame_columns():
    # #2113
    df = DataFrame({'x': [1, 2, 3], 'y': [3, 4, 5]})

    df.groupby(level=0, axis='columns').mean()
    df.groupby(level=0, axis='columns').mean()
    df.groupby(level=0, axis='columns').mean()
    df.groupby(level=0, axis='columns').mean()


def test_cython_agg_return_dict():
    # GH 16741
    df = DataFrame(
        {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
         'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
         'C': np.random.randn(8),
         'D': np.random.randn(8)})

    ts = df.groupby('A')['B'].agg(lambda x: x.value_counts().to_dict())
    expected = Series([{'two': 1, 'one': 1, 'three': 1},
                       {'two': 2, 'one': 2, 'three': 1}],
                      index=Index(['bar', 'foo'], name='A'),
                      name='B')
    tm.assert_series_equal(ts, expected)


def test_cython_fail_agg():
    dr = bdate_range('1/1/2000', periods=50)
    ts = Series(['A', 'B', 'C', 'D', 'E'] * 10, index=dr)

    grouped = ts.groupby(lambda x: x.month)
    summed = grouped.sum()
    expected = grouped.agg(np.sum)
    tm.assert_series_equal(summed, expected)


@pytest.mark.parametrize('op, targop', [
    ('mean', np.mean),
    ('median', np.median),
    ('var', np.var),
    ('add', np.sum),
    ('prod', np.prod),
    ('min', np.min),
    ('max', np.max),
    ('first', lambda x: x.iloc[0]),
    ('last', lambda x: x.iloc[-1]),
])
def test__cython_agg_general(op, targop):
    df = DataFrame(np.random.randn(1000))
    labels = np.random.randint(0, 50, size=1000).astype(float)

    result = df.groupby(labels)._cython_agg_general(op)
    expected = df.groupby(labels).agg(targop)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('op, targop', [
    ('mean', np.mean),
    ('median', lambda x: np.median(x) if len(x) > 0 else np.nan),
    ('var', lambda x: np.var(x, ddof=1)),
    ('min', np.min),
    ('max', np.max), ]
)
def test_cython_agg_empty_buckets(op, targop, observed):
    df = pd.DataFrame([11, 12, 13])
    grps = range(0, 55, 5)

    # calling _cython_agg_general directly, instead of via the user API
    # which sets different values for min_count, so do that here.
    g = df.groupby(pd.cut(df[0], grps), observed=observed)
    result = g._cython_agg_general(op)

    g = df.groupby(pd.cut(df[0], grps), observed=observed)
    expected = g.agg(lambda x: targop(x))
    tm.assert_frame_equal(result, expected)


def test_cython_agg_empty_buckets_nanops(observed):
    # GH-18869 can't call nanops on empty groups, so hardcode expected
    # for these
    df = pd.DataFrame([11, 12, 13], columns=['a'])
    grps = range(0, 25, 5)
    # add / sum
    result = df.groupby(pd.cut(df['a'], grps),
                        observed=observed)._cython_agg_general('add')
    intervals = pd.interval_range(0, 20, freq=5)
    expected = pd.DataFrame(
        {"a": [0, 0, 36, 0]},
        index=pd.CategoricalIndex(intervals, name='a', ordered=True))
    if observed:
        expected = expected[expected.a != 0]

    tm.assert_frame_equal(result, expected)

    # prod
    result = df.groupby(pd.cut(df['a'], grps),
                        observed=observed)._cython_agg_general('prod')
    expected = pd.DataFrame(
        {"a": [1, 1, 1716, 1]},
        index=pd.CategoricalIndex(intervals, name='a', ordered=True))
    if observed:
        expected = expected[expected.a != 1]

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('op', ['first', 'last', 'max', 'min'])
@pytest.mark.parametrize('data', [
    Timestamp('2016-10-14 21:00:44.557'),
    Timedelta('17088 days 21:00:44.557'), ])
def test_cython_with_timestamp_and_nat(op, data):
    # https://github.com/pandas-dev/pandas/issues/19526
    df = DataFrame({'a': [0, 1], 'b': [data, NaT]})
    index = Index([0, 1], name='a')

    # We will group by a and test the cython aggregations
    expected = DataFrame({'b': [data, NaT]}, index=index)

    result = df.groupby('a').aggregate(op)
    tm.assert_frame_equal(expected, result)
