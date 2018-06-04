import numpy as np
import pandas as pd
from pandas import DataFrame, MultiIndex, Index, Series, isna
from pandas.compat import lrange
from pandas.util.testing import (
    assert_frame_equal,
    assert_produces_warning,
    assert_series_equal)


def test_first_last_nth(df):
    # tests for first / last / nth
    grouped = df.groupby('A')
    first = grouped.first()
    expected = df.loc[[1, 0], ['B', 'C', 'D']]
    expected.index = Index(['bar', 'foo'], name='A')
    expected = expected.sort_index()
    assert_frame_equal(first, expected)

    nth = grouped.nth(0)
    assert_frame_equal(nth, expected)

    last = grouped.last()
    expected = df.loc[[5, 7], ['B', 'C', 'D']]
    expected.index = Index(['bar', 'foo'], name='A')
    assert_frame_equal(last, expected)

    nth = grouped.nth(-1)
    assert_frame_equal(nth, expected)

    nth = grouped.nth(1)
    expected = df.loc[[2, 3], ['B', 'C', 'D']].copy()
    expected.index = Index(['foo', 'bar'], name='A')
    expected = expected.sort_index()
    assert_frame_equal(nth, expected)

    # it works!
    grouped['B'].first()
    grouped['B'].last()
    grouped['B'].nth(0)

    df.loc[df['A'] == 'foo', 'B'] = np.nan
    assert isna(grouped['B'].first()['foo'])
    assert isna(grouped['B'].last()['foo'])
    assert isna(grouped['B'].nth(0)['foo'])

    # v0.14.0 whatsnew
    df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
    g = df.groupby('A')
    result = g.first()
    expected = df.iloc[[1, 2]].set_index('A')
    assert_frame_equal(result, expected)

    expected = df.iloc[[1, 2]].set_index('A')
    result = g.nth(0, dropna='any')
    assert_frame_equal(result, expected)


def test_first_last_nth_dtypes(df_mixed_floats):

    df = df_mixed_floats.copy()
    df['E'] = True
    df['F'] = 1

    # tests for first / last / nth
    grouped = df.groupby('A')
    first = grouped.first()
    expected = df.loc[[1, 0], ['B', 'C', 'D', 'E', 'F']]
    expected.index = Index(['bar', 'foo'], name='A')
    expected = expected.sort_index()
    assert_frame_equal(first, expected)

    last = grouped.last()
    expected = df.loc[[5, 7], ['B', 'C', 'D', 'E', 'F']]
    expected.index = Index(['bar', 'foo'], name='A')
    expected = expected.sort_index()
    assert_frame_equal(last, expected)

    nth = grouped.nth(1)
    expected = df.loc[[3, 2], ['B', 'C', 'D', 'E', 'F']]
    expected.index = Index(['bar', 'foo'], name='A')
    expected = expected.sort_index()
    assert_frame_equal(nth, expected)

    # GH 2763, first/last shifting dtypes
    idx = lrange(10)
    idx.append(9)
    s = Series(data=lrange(11), index=idx, name='IntCol')
    assert s.dtype == 'int64'
    f = s.groupby(level=0).first()
    assert f.dtype == 'int64'


def test_nth():
    df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
    g = df.groupby('A')

    assert_frame_equal(g.nth(0), df.iloc[[0, 2]].set_index('A'))
    assert_frame_equal(g.nth(1), df.iloc[[1]].set_index('A'))
    assert_frame_equal(g.nth(2), df.loc[[]].set_index('A'))
    assert_frame_equal(g.nth(-1), df.iloc[[1, 2]].set_index('A'))
    assert_frame_equal(g.nth(-2), df.iloc[[0]].set_index('A'))
    assert_frame_equal(g.nth(-3), df.loc[[]].set_index('A'))
    assert_series_equal(g.B.nth(0), df.set_index('A').B.iloc[[0, 2]])
    assert_series_equal(g.B.nth(1), df.set_index('A').B.iloc[[1]])
    assert_frame_equal(g[['B']].nth(0),
                       df.loc[[0, 2], ['A', 'B']].set_index('A'))

    exp = df.set_index('A')
    assert_frame_equal(g.nth(0, dropna='any'), exp.iloc[[1, 2]])
    assert_frame_equal(g.nth(-1, dropna='any'), exp.iloc[[1, 2]])

    exp['B'] = np.nan
    assert_frame_equal(g.nth(7, dropna='any'), exp.iloc[[1, 2]])
    assert_frame_equal(g.nth(2, dropna='any'), exp.iloc[[1, 2]])

    # out of bounds, regression from 0.13.1
    # GH 6621
    df = DataFrame({'color': {0: 'green',
                              1: 'green',
                              2: 'red',
                              3: 'red',
                              4: 'red'},
                    'food': {0: 'ham',
                             1: 'eggs',
                             2: 'eggs',
                             3: 'ham',
                             4: 'pork'},
                    'two': {0: 1.5456590000000001,
                            1: -0.070345000000000005,
                            2: -2.4004539999999999,
                            3: 0.46206000000000003,
                            4: 0.52350799999999997},
                    'one': {0: 0.56573799999999996,
                            1: -0.9742360000000001,
                            2: 1.033801,
                            3: -0.78543499999999999,
                            4: 0.70422799999999997}}).set_index(['color',
                                                                 'food'])

    result = df.groupby(level=0, as_index=False).nth(2)
    expected = df.iloc[[-1]]
    assert_frame_equal(result, expected)

    result = df.groupby(level=0, as_index=False).nth(3)
    expected = df.loc[[]]
    assert_frame_equal(result, expected)

    # GH 7559
    # from the vbench
    df = DataFrame(np.random.randint(1, 10, (100, 2)), dtype='int64')
    s = df[1]
    g = df[0]
    expected = s.groupby(g).first()
    expected2 = s.groupby(g).apply(lambda x: x.iloc[0])
    assert_series_equal(expected2, expected, check_names=False)
    assert expected.name == 1
    assert expected2.name == 1

    # validate first
    v = s[g == 1].iloc[0]
    assert expected.iloc[0] == v
    assert expected2.iloc[0] == v

    # this is NOT the same as .first (as sorted is default!)
    # as it keeps the order in the series (and not the group order)
    # related GH 7287
    expected = s.groupby(g, sort=False).first()
    result = s.groupby(g, sort=False).nth(0, dropna='all')
    assert_series_equal(result, expected)

    # doc example
    df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
    g = df.groupby('A')
    # PR 17493, related to issue 11038
    # test Series.nth with True for dropna produces FutureWarning
    with assert_produces_warning(FutureWarning):
        result = g.B.nth(0, dropna=True)
    expected = g.B.first()
    assert_series_equal(result, expected)

    # test multiple nth values
    df = DataFrame([[1, np.nan], [1, 3], [1, 4], [5, 6], [5, 7]],
                   columns=['A', 'B'])
    g = df.groupby('A')

    assert_frame_equal(g.nth(0), df.iloc[[0, 3]].set_index('A'))
    assert_frame_equal(g.nth([0]), df.iloc[[0, 3]].set_index('A'))
    assert_frame_equal(g.nth([0, 1]), df.iloc[[0, 1, 3, 4]].set_index('A'))
    assert_frame_equal(
        g.nth([0, -1]), df.iloc[[0, 2, 3, 4]].set_index('A'))
    assert_frame_equal(
        g.nth([0, 1, 2]), df.iloc[[0, 1, 2, 3, 4]].set_index('A'))
    assert_frame_equal(
        g.nth([0, 1, -1]), df.iloc[[0, 1, 2, 3, 4]].set_index('A'))
    assert_frame_equal(g.nth([2]), df.iloc[[2]].set_index('A'))
    assert_frame_equal(g.nth([3, 4]), df.loc[[]].set_index('A'))

    business_dates = pd.date_range(start='4/1/2014', end='6/30/2014',
                                   freq='B')
    df = DataFrame(1, index=business_dates, columns=['a', 'b'])
    # get the first, fourth and last two business days for each month
    key = [df.index.year, df.index.month]
    result = df.groupby(key, as_index=False).nth([0, 3, -2, -1])
    expected_dates = pd.to_datetime(
        ['2014/4/1', '2014/4/4', '2014/4/29', '2014/4/30', '2014/5/1',
         '2014/5/6', '2014/5/29', '2014/5/30', '2014/6/2', '2014/6/5',
         '2014/6/27', '2014/6/30'])
    expected = DataFrame(1, columns=['a', 'b'], index=expected_dates)
    assert_frame_equal(result, expected)


def test_nth_multi_index(three_group):
    # PR 9090, related to issue 8979
    # test nth on MultiIndex, should match .first()
    grouped = three_group.groupby(['A', 'B'])
    result = grouped.nth(0)
    expected = grouped.first()
    assert_frame_equal(result, expected)


def test_nth_multi_index_as_expected():
    # PR 9090, related to issue 8979
    # test nth on MultiIndex
    three_group = DataFrame(
        {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
               'foo', 'foo', 'foo'],
         'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
               'two', 'two', 'one'],
         'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
               'dull', 'shiny', 'shiny', 'shiny']})
    grouped = three_group.groupby(['A', 'B'])
    result = grouped.nth(0)
    expected = DataFrame(
        {'C': ['dull', 'dull', 'dull', 'dull']},
        index=MultiIndex.from_arrays([['bar', 'bar', 'foo', 'foo'],
                                      ['one', 'two', 'one', 'two']],
                                     names=['A', 'B']))
    assert_frame_equal(result, expected)


def test_groupby_head_tail():
    df = DataFrame([[1, 2], [1, 4], [5, 6]], columns=['A', 'B'])
    g_as = df.groupby('A', as_index=True)
    g_not_as = df.groupby('A', as_index=False)

    # as_index= False, much easier
    assert_frame_equal(df.loc[[0, 2]], g_not_as.head(1))
    assert_frame_equal(df.loc[[1, 2]], g_not_as.tail(1))

    empty_not_as = DataFrame(columns=df.columns,
                             index=pd.Index([], dtype=df.index.dtype))
    empty_not_as['A'] = empty_not_as['A'].astype(df.A.dtype)
    empty_not_as['B'] = empty_not_as['B'].astype(df.B.dtype)
    assert_frame_equal(empty_not_as, g_not_as.head(0))
    assert_frame_equal(empty_not_as, g_not_as.tail(0))
    assert_frame_equal(empty_not_as, g_not_as.head(-1))
    assert_frame_equal(empty_not_as, g_not_as.tail(-1))

    assert_frame_equal(df, g_not_as.head(7))  # contains all
    assert_frame_equal(df, g_not_as.tail(7))

    # as_index=True, (used to be different)
    df_as = df

    assert_frame_equal(df_as.loc[[0, 2]], g_as.head(1))
    assert_frame_equal(df_as.loc[[1, 2]], g_as.tail(1))

    empty_as = DataFrame(index=df_as.index[:0], columns=df.columns)
    empty_as['A'] = empty_not_as['A'].astype(df.A.dtype)
    empty_as['B'] = empty_not_as['B'].astype(df.B.dtype)
    assert_frame_equal(empty_as, g_as.head(0))
    assert_frame_equal(empty_as, g_as.tail(0))
    assert_frame_equal(empty_as, g_as.head(-1))
    assert_frame_equal(empty_as, g_as.tail(-1))

    assert_frame_equal(df_as, g_as.head(7))  # contains all
    assert_frame_equal(df_as, g_as.tail(7))

    # test with selection
    assert_frame_equal(g_as[[]].head(1), df_as.loc[[0, 2], []])
    assert_frame_equal(g_as[['A']].head(1), df_as.loc[[0, 2], ['A']])
    assert_frame_equal(g_as[['B']].head(1), df_as.loc[[0, 2], ['B']])
    assert_frame_equal(g_as[['A', 'B']].head(1), df_as.loc[[0, 2]])

    assert_frame_equal(g_not_as[[]].head(1), df_as.loc[[0, 2], []])
    assert_frame_equal(g_not_as[['A']].head(1), df_as.loc[[0, 2], ['A']])
    assert_frame_equal(g_not_as[['B']].head(1), df_as.loc[[0, 2], ['B']])
    assert_frame_equal(g_not_as[['A', 'B']].head(1), df_as.loc[[0, 2]])


def test_group_selection_cache():
    # GH 12839 nth, head, and tail should return same result consistently
    df = DataFrame([[1, 2], [1, 4], [5, 6]], columns=['A', 'B'])
    expected = df.iloc[[0, 2]].set_index('A')

    g = df.groupby('A')
    result1 = g.head(n=2)
    result2 = g.nth(0)
    assert_frame_equal(result1, df)
    assert_frame_equal(result2, expected)

    g = df.groupby('A')
    result1 = g.tail(n=2)
    result2 = g.nth(0)
    assert_frame_equal(result1, df)
    assert_frame_equal(result2, expected)

    g = df.groupby('A')
    result1 = g.nth(0)
    result2 = g.head(n=2)
    assert_frame_equal(result1, expected)
    assert_frame_equal(result2, df)

    g = df.groupby('A')
    result1 = g.nth(0)
    result2 = g.tail(n=2)
    assert_frame_equal(result1, expected)
    assert_frame_equal(result2, df)


def test_nth_empty():
    # GH 16064
    df = DataFrame(index=[0], columns=['a', 'b', 'c'])
    result = df.groupby('a').nth(10)
    expected = DataFrame(index=Index([], name='a'), columns=['b', 'c'])
    assert_frame_equal(result, expected)

    result = df.groupby(['a', 'b']).nth(10)
    expected = DataFrame(index=MultiIndex([[], []], [[], []],
                                          names=['a', 'b']),
                         columns=['c'])
    assert_frame_equal(result, expected)
