import pytest
import numpy as np
import pandas as pd
from datetime import datetime
from pandas.util import testing as tm
from pandas import DataFrame, MultiIndex, compat, Series, bdate_range, Index


def test_apply_issues():
        # GH 5788

    s = """2011.05.16,00:00,1.40893
2011.05.16,01:00,1.40760
2011.05.16,02:00,1.40750
2011.05.16,03:00,1.40649
2011.05.17,02:00,1.40893
2011.05.17,03:00,1.40760
2011.05.17,04:00,1.40750
2011.05.17,05:00,1.40649
2011.05.18,02:00,1.40893
2011.05.18,03:00,1.40760
2011.05.18,04:00,1.40750
2011.05.18,05:00,1.40649"""

    df = pd.read_csv(
        compat.StringIO(s), header=None, names=['date', 'time', 'value'],
        parse_dates=[['date', 'time']])
    df = df.set_index('date_time')

    expected = df.groupby(df.index.date).idxmax()
    result = df.groupby(df.index.date).apply(lambda x: x.idxmax())
    tm.assert_frame_equal(result, expected)

    # GH 5789
    # don't auto coerce dates
    df = pd.read_csv(
        compat.StringIO(s), header=None, names=['date', 'time', 'value'])
    exp_idx = pd.Index(
        ['2011.05.16', '2011.05.17', '2011.05.18'
         ], dtype=object, name='date')
    expected = Series(['00:00', '02:00', '02:00'], index=exp_idx)
    result = df.groupby('date').apply(
        lambda x: x['time'][x['value'].idxmax()])
    tm.assert_series_equal(result, expected)


def test_apply_trivial():
    # GH 20066
    # trivial apply: ignore input and return a constant dataframe.
    df = pd.DataFrame({'key': ['a', 'a', 'b', 'b', 'a'],
                       'data': [1.0, 2.0, 3.0, 4.0, 5.0]},
                      columns=['key', 'data'])
    expected = pd.concat([df.iloc[1:], df.iloc[1:]],
                         axis=1, keys=['float64', 'object'])
    result = df.groupby([str(x) for x in df.dtypes],
                        axis=1).apply(lambda x: df.iloc[1:])

    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason=("GH 20066; function passed into apply "
                           "returns a DataFrame with the same index "
                           "as the one to create GroupBy object."))
def test_apply_trivial_fail():
    # GH 20066
    # trivial apply fails if the constant dataframe has the same index
    # with the one used to create GroupBy object.
    df = pd.DataFrame({'key': ['a', 'a', 'b', 'b', 'a'],
                       'data': [1.0, 2.0, 3.0, 4.0, 5.0]},
                      columns=['key', 'data'])
    expected = pd.concat([df, df],
                         axis=1, keys=['float64', 'object'])
    result = df.groupby([str(x) for x in df.dtypes],
                        axis=1).apply(lambda x: df)

    tm.assert_frame_equal(result, expected)


def test_fast_apply():
    # make sure that fast apply is correctly called
    # rather than raising any kind of error
    # otherwise the python path will be callsed
    # which slows things down
    N = 1000
    labels = np.random.randint(0, 2000, size=N)
    labels2 = np.random.randint(0, 3, size=N)
    df = DataFrame({'key': labels,
                    'key2': labels2,
                    'value1': np.random.randn(N),
                    'value2': ['foo', 'bar', 'baz', 'qux'] * (N // 4)})

    def f(g):
        return 1

    g = df.groupby(['key', 'key2'])

    grouper = g.grouper

    splitter = grouper._get_splitter(g._selected_obj, axis=g.axis)
    group_keys = grouper._get_group_keys()

    values, mutated = splitter.fast_apply(f, group_keys)
    assert not mutated


def test_apply_with_mixed_dtype():
    # GH3480, apply with mixed dtype on axis=1 breaks in 0.11
    df = DataFrame({'foo1': np.random.randn(6),
                    'foo2': ['one', 'two', 'two', 'three', 'one', 'two']})
    result = df.apply(lambda x: x, axis=1)
    tm.assert_series_equal(df.get_dtype_counts(), result.get_dtype_counts())

    # GH 3610 incorrect dtype conversion with as_index=False
    df = DataFrame({"c1": [1, 2, 6, 6, 8]})
    df["c2"] = df.c1 / 2.0
    result1 = df.groupby("c2").mean().reset_index().c2
    result2 = df.groupby("c2", as_index=False).mean().c2
    tm.assert_series_equal(result1, result2)


def test_groupby_as_index_apply(df):
    # GH #4648 and #3417
    df = DataFrame({'item_id': ['b', 'b', 'a', 'c', 'a', 'b'],
                    'user_id': [1, 2, 1, 1, 3, 1],
                    'time': range(6)})

    g_as = df.groupby('user_id', as_index=True)
    g_not_as = df.groupby('user_id', as_index=False)

    res_as = g_as.head(2).index
    res_not_as = g_not_as.head(2).index
    exp = Index([0, 1, 2, 4])
    tm.assert_index_equal(res_as, exp)
    tm.assert_index_equal(res_not_as, exp)

    res_as_apply = g_as.apply(lambda x: x.head(2)).index
    res_not_as_apply = g_not_as.apply(lambda x: x.head(2)).index

    # apply doesn't maintain the original ordering
    # changed in GH5610 as the as_index=False returns a MI here
    exp_not_as_apply = MultiIndex.from_tuples([(0, 0), (0, 2), (1, 1), (
        2, 4)])
    tp = [(1, 0), (1, 2), (2, 1), (3, 4)]
    exp_as_apply = MultiIndex.from_tuples(tp, names=['user_id', None])

    tm.assert_index_equal(res_as_apply, exp_as_apply)
    tm.assert_index_equal(res_not_as_apply, exp_not_as_apply)

    ind = Index(list('abcde'))
    df = DataFrame([[1, 2], [2, 3], [1, 4], [1, 5], [2, 6]], index=ind)
    res = df.groupby(0, as_index=False).apply(lambda x: x).index
    tm.assert_index_equal(res, ind)


def test_apply_concat_preserve_names(three_group):
    grouped = three_group.groupby(['A', 'B'])

    def desc(group):
        result = group.describe()
        result.index.name = 'stat'
        return result

    def desc2(group):
        result = group.describe()
        result.index.name = 'stat'
        result = result[:len(group)]
        # weirdo
        return result

    def desc3(group):
        result = group.describe()

        # names are different
        result.index.name = 'stat_%d' % len(group)

        result = result[:len(group)]
        # weirdo
        return result

    result = grouped.apply(desc)
    assert result.index.names == ('A', 'B', 'stat')

    result2 = grouped.apply(desc2)
    assert result2.index.names == ('A', 'B', 'stat')

    result3 = grouped.apply(desc3)
    assert result3.index.names == ('A', 'B', None)


def test_apply_series_to_frame():
    def f(piece):
        with np.errstate(invalid='ignore'):
            logged = np.log(piece)
        return DataFrame({'value': piece,
                          'demeaned': piece - piece.mean(),
                          'logged': logged})

    dr = bdate_range('1/1/2000', periods=100)
    ts = Series(np.random.randn(100), index=dr)

    grouped = ts.groupby(lambda x: x.month)
    result = grouped.apply(f)

    assert isinstance(result, DataFrame)
    tm.assert_index_equal(result.index, ts.index)


def test_apply_series_yield_constant(df):
    result = df.groupby(['A', 'B'])['C'].apply(len)
    assert result.index.names[:2] == ('A', 'B')


def test_apply_frame_yield_constant(df):
    # GH13568
    result = df.groupby(['A', 'B']).apply(len)
    assert isinstance(result, Series)
    assert result.name is None

    result = df.groupby(['A', 'B'])[['C', 'D']].apply(len)
    assert isinstance(result, Series)
    assert result.name is None


def test_apply_frame_to_series(df):
    grouped = df.groupby(['A', 'B'])
    result = grouped.apply(len)
    expected = grouped.count()['C']
    tm.assert_index_equal(result.index, expected.index)
    tm.assert_numpy_array_equal(result.values, expected.values)


def test_apply_frame_concat_series():
    def trans(group):
        return group.groupby('B')['C'].sum().sort_values()[:2]

    def trans2(group):
        grouped = group.groupby(df.reindex(group.index)['B'])
        return grouped.sum().sort_values()[:2]

    df = DataFrame({'A': np.random.randint(0, 5, 1000),
                    'B': np.random.randint(0, 5, 1000),
                    'C': np.random.randn(1000)})

    result = df.groupby('A').apply(trans)
    exp = df.groupby('A')['C'].apply(trans2)
    tm.assert_series_equal(result, exp, check_names=False)
    assert result.name == 'C'


def test_apply_transform(ts):
    grouped = ts.groupby(lambda x: x.month)
    result = grouped.apply(lambda x: x * 2)
    expected = grouped.transform(lambda x: x * 2)
    tm.assert_series_equal(result, expected)


def test_apply_multikey_corner(tsframe):
    grouped = tsframe.groupby([lambda x: x.year, lambda x: x.month])

    def f(group):
        return group.sort_values('A')[-5:]

    result = grouped.apply(f)
    for key, group in grouped:
        tm.assert_frame_equal(result.loc[key], f(group))


def test_apply_chunk_view():
    # Low level tinkering could be unsafe, make sure not
    df = DataFrame({'key': [1, 1, 1, 2, 2, 2, 3, 3, 3],
                    'value': compat.lrange(9)})

    # return view
    f = lambda x: x[:2]

    result = df.groupby('key', group_keys=False).apply(f)
    expected = df.take([0, 1, 3, 4, 6, 7])
    tm.assert_frame_equal(result, expected)


def test_apply_no_name_column_conflict():
    df = DataFrame({'name': [1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
                    'name2': [0, 0, 0, 1, 1, 1, 0, 0, 1, 1],
                    'value': compat.lrange(10)[::-1]})

    # it works! #2605
    grouped = df.groupby(['name', 'name2'])
    grouped.apply(lambda x: x.sort_values('value', inplace=True))


def test_apply_typecast_fail():
    df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                    'c': np.tile(
                        ['a', 'b', 'c'], 2),
                    'v': np.arange(1., 7.)})

    def f(group):
        v = group['v']
        group['v2'] = (v - v.min()) / (v.max() - v.min())
        return group

    result = df.groupby('d').apply(f)

    expected = df.copy()
    expected['v2'] = np.tile([0., 0.5, 1], 2)

    tm.assert_frame_equal(result, expected)


def test_apply_multiindex_fail():
    index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [1, 2, 3, 1, 2, 3]
                                    ])
    df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                    'c': np.tile(['a', 'b', 'c'], 2),
                    'v': np.arange(1., 7.)}, index=index)

    def f(group):
        v = group['v']
        group['v2'] = (v - v.min()) / (v.max() - v.min())
        return group

    result = df.groupby('d').apply(f)

    expected = df.copy()
    expected['v2'] = np.tile([0., 0.5, 1], 2)

    tm.assert_frame_equal(result, expected)


def test_apply_corner(tsframe):
    result = tsframe.groupby(lambda x: x.year).apply(lambda x: x * 2)
    expected = tsframe * 2
    tm.assert_frame_equal(result, expected)


def test_apply_without_copy():
    # GH 5545
    # returning a non-copy in an applied function fails

    data = DataFrame({'id_field': [100, 100, 200, 300],
                      'category': ['a', 'b', 'c', 'c'],
                      'value': [1, 2, 3, 4]})

    def filt1(x):
        if x.shape[0] == 1:
            return x.copy()
        else:
            return x[x.category == 'c']

    def filt2(x):
        if x.shape[0] == 1:
            return x
        else:
            return x[x.category == 'c']

    expected = data.groupby('id_field').apply(filt1)
    result = data.groupby('id_field').apply(filt2)
    tm.assert_frame_equal(result, expected)


def test_apply_corner_cases():
    # #535, can't use sliding iterator

    N = 1000
    labels = np.random.randint(0, 100, size=N)
    df = DataFrame({'key': labels,
                    'value1': np.random.randn(N),
                    'value2': ['foo', 'bar', 'baz', 'qux'] * (N // 4)})

    grouped = df.groupby('key')

    def f(g):
        g['value3'] = g['value1'] * 2
        return g

    result = grouped.apply(f)
    assert 'value3' in result


def test_apply_numeric_coercion_when_datetime():
    # In the past, group-by/apply operations have been over-eager
    # in converting dtypes to numeric, in the presence of datetime
    # columns.  Various GH issues were filed, the reproductions
    # for which are here.

    # GH 15670
    df = pd.DataFrame({'Number': [1, 2],
                       'Date': ["2017-03-02"] * 2,
                       'Str': ["foo", "inf"]})
    expected = df.groupby(['Number']).apply(lambda x: x.iloc[0])
    df.Date = pd.to_datetime(df.Date)
    result = df.groupby(['Number']).apply(lambda x: x.iloc[0])
    tm.assert_series_equal(result['Str'], expected['Str'])

    # GH 15421
    df = pd.DataFrame({'A': [10, 20, 30],
                       'B': ['foo', '3', '4'],
                       'T': [pd.Timestamp("12:31:22")] * 3})

    def get_B(g):
        return g.iloc[0][['B']]
    result = df.groupby('A').apply(get_B)['B']
    expected = df.B
    expected.index = df.A
    tm.assert_series_equal(result, expected)

    # GH 14423
    def predictions(tool):
        out = pd.Series(index=['p1', 'p2', 'useTime'], dtype=object)
        if 'step1' in list(tool.State):
            out['p1'] = str(tool[tool.State == 'step1'].Machine.values[0])
        if 'step2' in list(tool.State):
            out['p2'] = str(tool[tool.State == 'step2'].Machine.values[0])
            out['useTime'] = str(
                tool[tool.State == 'step2'].oTime.values[0])
        return out
    df1 = pd.DataFrame({'Key': ['B', 'B', 'A', 'A'],
                        'State': ['step1', 'step2', 'step1', 'step2'],
                        'oTime': ['', '2016-09-19 05:24:33',
                                  '', '2016-09-19 23:59:04'],
                        'Machine': ['23', '36L', '36R', '36R']})
    df2 = df1.copy()
    df2.oTime = pd.to_datetime(df2.oTime)
    expected = df1.groupby('Key').apply(predictions).p1
    result = df2.groupby('Key').apply(predictions).p1
    tm.assert_series_equal(expected, result)


def test_time_field_bug():
    # Test a fix for the following error related to GH issue 11324 When
    # non-key fields in a group-by dataframe contained time-based fields
    # that were not returned by the apply function, an exception would be
    # raised.

    df = pd.DataFrame({'a': 1, 'b': [datetime.now() for nn in range(10)]})

    def func_with_no_date(batch):
        return pd.Series({'c': 2})

    def func_with_date(batch):
        return pd.Series({'b': datetime(2015, 1, 1), 'c': 2})

    dfg_no_conversion = df.groupby(by=['a']).apply(func_with_no_date)
    dfg_no_conversion_expected = pd.DataFrame({'c': 2}, index=[1])
    dfg_no_conversion_expected.index.name = 'a'

    dfg_conversion = df.groupby(by=['a']).apply(func_with_date)
    dfg_conversion_expected = pd.DataFrame(
        {'b': datetime(2015, 1, 1),
         'c': 2}, index=[1])
    dfg_conversion_expected.index.name = 'a'

    tm.assert_frame_equal(dfg_no_conversion, dfg_no_conversion_expected)
    tm.assert_frame_equal(dfg_conversion, dfg_conversion_expected)


def test_gb_apply_list_of_unequal_len_arrays():

    # GH1738
    df = DataFrame({'group1': ['a', 'a', 'a', 'b', 'b', 'b', 'a', 'a', 'a',
                               'b', 'b', 'b'],
                    'group2': ['c', 'c', 'd', 'd', 'd', 'e', 'c', 'c', 'd',
                               'd', 'd', 'e'],
                    'weight': [1.1, 2, 3, 4, 5, 6, 2, 4, 6, 8, 1, 2],
                    'value': [7.1, 8, 9, 10, 11, 12, 8, 7, 6, 5, 4, 3]})
    df = df.set_index(['group1', 'group2'])
    df_grouped = df.groupby(level=['group1', 'group2'], sort=True)

    def noddy(value, weight):
        out = np.array(value * weight).repeat(3)
        return out

    # the kernel function returns arrays of unequal length
    # pandas sniffs the first one, sees it's an array and not
    # a list, and assumed the rest are of equal length
    # and so tries a vstack

    # don't die
    df_grouped.apply(lambda x: noddy(x.value, x.weight))


def test_groupby_apply_all_none():
    # Tests to make sure no errors if apply function returns all None
    # values. Issue 9684.
    test_df = DataFrame({'groups': [0, 0, 1, 1],
                         'random_vars': [8, 7, 4, 5]})

    def test_func(x):
        pass

    result = test_df.groupby('groups').apply(test_func)
    expected = DataFrame()
    tm.assert_frame_equal(result, expected)


def test_groupby_apply_none_first():
    # GH 12824. Tests if apply returns None first.
    test_df1 = DataFrame({'groups': [1, 1, 1, 2], 'vars': [0, 1, 2, 3]})
    test_df2 = DataFrame({'groups': [1, 2, 2, 2], 'vars': [0, 1, 2, 3]})

    def test_func(x):
        if x.shape[0] < 2:
            return None
        return x.iloc[[0, -1]]

    result1 = test_df1.groupby('groups').apply(test_func)
    result2 = test_df2.groupby('groups').apply(test_func)
    index1 = MultiIndex.from_arrays([[1, 1], [0, 2]],
                                    names=['groups', None])
    index2 = MultiIndex.from_arrays([[2, 2], [1, 3]],
                                    names=['groups', None])
    expected1 = DataFrame({'groups': [1, 1], 'vars': [0, 2]},
                          index=index1)
    expected2 = DataFrame({'groups': [2, 2], 'vars': [1, 3]},
                          index=index2)
    tm.assert_frame_equal(result1, expected1)
    tm.assert_frame_equal(result2, expected2)


def test_apply_with_mixed_types():
    # gh-20949
    df = pd.DataFrame({'A': 'a a b'.split(), 'B': [1, 2, 3], 'C': [4, 6, 5]})
    g = df.groupby('A')

    result = g.transform(lambda x: x / x.sum())
    expected = pd.DataFrame({'B': [1 / 3., 2 / 3., 1], 'C': [0.4, 0.6, 1.0]})
    tm.assert_frame_equal(result, expected)

    result = g.apply(lambda x: x / x.sum())
    tm.assert_frame_equal(result, expected)
