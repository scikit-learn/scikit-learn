import pytest
import numpy as np
import pandas as pd
from pandas import DataFrame, concat
from pandas.util import testing as tm


def test_rank_apply():
    lev1 = tm.rands_array(10, 100)
    lev2 = tm.rands_array(10, 130)
    lab1 = np.random.randint(0, 100, size=500)
    lab2 = np.random.randint(0, 130, size=500)

    df = DataFrame({'value': np.random.randn(500),
                    'key1': lev1.take(lab1),
                    'key2': lev2.take(lab2)})

    result = df.groupby(['key1', 'key2']).value.rank()

    expected = []
    for key, piece in df.groupby(['key1', 'key2']):
        expected.append(piece.value.rank())
    expected = concat(expected, axis=0)
    expected = expected.reindex(result.index)
    tm.assert_series_equal(result, expected)

    result = df.groupby(['key1', 'key2']).value.rank(pct=True)

    expected = []
    for key, piece in df.groupby(['key1', 'key2']):
        expected.append(piece.value.rank(pct=True))
    expected = concat(expected, axis=0)
    expected = expected.reindex(result.index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("grps", [
    ['qux'], ['qux', 'quux']])
@pytest.mark.parametrize("vals", [
    [2, 2, 8, 2, 6],
    [pd.Timestamp('2018-01-02'), pd.Timestamp('2018-01-02'),
     pd.Timestamp('2018-01-08'), pd.Timestamp('2018-01-02'),
     pd.Timestamp('2018-01-06')]])
@pytest.mark.parametrize("ties_method,ascending,pct,exp", [
    ('average', True, False, [2., 2., 5., 2., 4.]),
    ('average', True, True, [0.4, 0.4, 1.0, 0.4, 0.8]),
    ('average', False, False, [4., 4., 1., 4., 2.]),
    ('average', False, True, [.8, .8, .2, .8, .4]),
    ('min', True, False, [1., 1., 5., 1., 4.]),
    ('min', True, True, [0.2, 0.2, 1.0, 0.2, 0.8]),
    ('min', False, False, [3., 3., 1., 3., 2.]),
    ('min', False, True, [.6, .6, .2, .6, .4]),
    ('max', True, False, [3., 3., 5., 3., 4.]),
    ('max', True, True, [0.6, 0.6, 1.0, 0.6, 0.8]),
    ('max', False, False, [5., 5., 1., 5., 2.]),
    ('max', False, True, [1., 1., .2, 1., .4]),
    ('first', True, False, [1., 2., 5., 3., 4.]),
    ('first', True, True, [0.2, 0.4, 1.0, 0.6, 0.8]),
    ('first', False, False, [3., 4., 1., 5., 2.]),
    ('first', False, True, [.6, .8, .2, 1., .4]),
    ('dense', True, False, [1., 1., 3., 1., 2.]),
    ('dense', True, True, [0.2, 0.2, 0.6, 0.2, 0.4]),
    ('dense', False, False, [3., 3., 1., 3., 2.]),
    ('dense', False, True, [.6, .6, .2, .6, .4]),
])
def test_rank_args(grps, vals, ties_method, ascending, pct, exp):
    key = np.repeat(grps, len(vals))
    vals = vals * len(grps)
    df = DataFrame({'key': key, 'val': vals})
    result = df.groupby('key').rank(method=ties_method,
                                    ascending=ascending, pct=pct)

    exp_df = DataFrame(exp * len(grps), columns=['val'])
    tm.assert_frame_equal(result, exp_df)


@pytest.mark.parametrize("grps", [
    ['qux'], ['qux', 'quux']])
@pytest.mark.parametrize("vals", [
    [-np.inf, -np.inf, np.nan, 1., np.nan, np.inf, np.inf],
])
@pytest.mark.parametrize("ties_method,ascending,na_option,exp", [
    ('average', True, 'keep', [1.5, 1.5, np.nan, 3, np.nan, 4.5, 4.5]),
    ('average', True, 'top', [3.5, 3.5, 1.5, 5., 1.5, 6.5, 6.5]),
    ('average', True, 'bottom', [1.5, 1.5, 6.5, 3., 6.5, 4.5, 4.5]),
    ('average', False, 'keep', [4.5, 4.5, np.nan, 3, np.nan, 1.5, 1.5]),
    ('average', False, 'top', [6.5, 6.5, 1.5, 5., 1.5, 3.5, 3.5]),
    ('average', False, 'bottom', [4.5, 4.5, 6.5, 3., 6.5, 1.5, 1.5]),
    ('min', True, 'keep', [1., 1., np.nan, 3., np.nan, 4., 4.]),
    ('min', True, 'top', [3., 3., 1., 5., 1., 6., 6.]),
    ('min', True, 'bottom', [1., 1., 6., 3., 6., 4., 4.]),
    ('min', False, 'keep', [4., 4., np.nan, 3., np.nan, 1., 1.]),
    ('min', False, 'top', [6., 6., 1., 5., 1., 3., 3.]),
    ('min', False, 'bottom', [4., 4., 6., 3., 6., 1., 1.]),
    ('max', True, 'keep', [2., 2., np.nan, 3., np.nan, 5., 5.]),
    ('max', True, 'top', [4., 4., 2., 5., 2., 7., 7.]),
    ('max', True, 'bottom', [2., 2., 7., 3., 7., 5., 5.]),
    ('max', False, 'keep', [5., 5., np.nan, 3., np.nan, 2., 2.]),
    ('max', False, 'top', [7., 7., 2., 5., 2., 4., 4.]),
    ('max', False, 'bottom', [5., 5., 7., 3., 7., 2., 2.]),
    ('first', True, 'keep', [1., 2., np.nan, 3., np.nan, 4., 5.]),
    ('first', True, 'top', [3., 4., 1., 5., 2., 6., 7.]),
    ('first', True, 'bottom', [1., 2., 6., 3., 7., 4., 5.]),
    ('first', False, 'keep', [4., 5., np.nan, 3., np.nan, 1., 2.]),
    ('first', False, 'top', [6., 7., 1., 5., 2., 3., 4.]),
    ('first', False, 'bottom', [4., 5., 6., 3., 7., 1., 2.]),
    ('dense', True, 'keep', [1., 1., np.nan, 2., np.nan, 3., 3.]),
    ('dense', True, 'top', [2., 2., 1., 3., 1., 4., 4.]),
    ('dense', True, 'bottom', [1., 1., 4., 2., 4., 3., 3.]),
    ('dense', False, 'keep', [3., 3., np.nan, 2., np.nan, 1., 1.]),
    ('dense', False, 'top', [4., 4., 1., 3., 1., 2., 2.]),
    ('dense', False, 'bottom', [3., 3., 4., 2., 4., 1., 1.])
])
def test_infs_n_nans(grps, vals, ties_method, ascending, na_option, exp):
    # GH 20561
    key = np.repeat(grps, len(vals))
    vals = vals * len(grps)
    df = DataFrame({'key': key, 'val': vals})
    result = df.groupby('key').rank(method=ties_method,
                                    ascending=ascending,
                                    na_option=na_option)
    exp_df = DataFrame(exp * len(grps), columns=['val'])
    tm.assert_frame_equal(result, exp_df)


@pytest.mark.parametrize("grps", [
    ['qux'], ['qux', 'quux']])
@pytest.mark.parametrize("vals", [
    [2, 2, np.nan, 8, 2, 6, np.nan, np.nan],  # floats
    [pd.Timestamp('2018-01-02'), pd.Timestamp('2018-01-02'), np.nan,
     pd.Timestamp('2018-01-08'), pd.Timestamp('2018-01-02'),
     pd.Timestamp('2018-01-06'), np.nan, np.nan]
])
@pytest.mark.parametrize("ties_method,ascending,na_option,pct,exp", [
    ('average', True, 'keep', False,
        [2., 2., np.nan, 5., 2., 4., np.nan, np.nan]),
    ('average', True, 'keep', True,
        [0.4, 0.4, np.nan, 1.0, 0.4, 0.8, np.nan, np.nan]),
    ('average', False, 'keep', False,
        [4., 4., np.nan, 1., 4., 2., np.nan, np.nan]),
    ('average', False, 'keep', True,
        [.8, 0.8, np.nan, 0.2, 0.8, 0.4, np.nan, np.nan]),
    ('min', True, 'keep', False,
        [1., 1., np.nan, 5., 1., 4., np.nan, np.nan]),
    ('min', True, 'keep', True,
        [0.2, 0.2, np.nan, 1.0, 0.2, 0.8, np.nan, np.nan]),
    ('min', False, 'keep', False,
        [3., 3., np.nan, 1., 3., 2., np.nan, np.nan]),
    ('min', False, 'keep', True,
        [.6, 0.6, np.nan, 0.2, 0.6, 0.4, np.nan, np.nan]),
    ('max', True, 'keep', False,
        [3., 3., np.nan, 5., 3., 4., np.nan, np.nan]),
    ('max', True, 'keep', True,
        [0.6, 0.6, np.nan, 1.0, 0.6, 0.8, np.nan, np.nan]),
    ('max', False, 'keep', False,
        [5., 5., np.nan, 1., 5., 2., np.nan, np.nan]),
    ('max', False, 'keep', True,
        [1., 1., np.nan, 0.2, 1., 0.4, np.nan, np.nan]),
    ('first', True, 'keep', False,
        [1., 2., np.nan, 5., 3., 4., np.nan, np.nan]),
    ('first', True, 'keep', True,
        [0.2, 0.4, np.nan, 1.0, 0.6, 0.8, np.nan, np.nan]),
    ('first', False, 'keep', False,
        [3., 4., np.nan, 1., 5., 2., np.nan, np.nan]),
    ('first', False, 'keep', True,
        [.6, 0.8, np.nan, 0.2, 1., 0.4, np.nan, np.nan]),
    ('dense', True, 'keep', False,
        [1., 1., np.nan, 3., 1., 2., np.nan, np.nan]),
    ('dense', True, 'keep', True,
        [0.2, 0.2, np.nan, 0.6, 0.2, 0.4, np.nan, np.nan]),
    ('dense', False, 'keep', False,
        [3., 3., np.nan, 1., 3., 2., np.nan, np.nan]),
    ('dense', False, 'keep', True,
        [.6, 0.6, np.nan, 0.2, 0.6, 0.4, np.nan, np.nan]),
    ('average', True, 'no_na', False, [2., 2., 7., 5., 2., 4., 7., 7.]),
    ('average', True, 'no_na', True,
        [0.25, 0.25, 0.875, 0.625, 0.25, 0.5, 0.875, 0.875]),
    ('average', False, 'no_na', False, [4., 4., 7., 1., 4., 2., 7., 7.]),
    ('average', False, 'no_na', True,
        [0.5, 0.5, 0.875, 0.125, 0.5, 0.25, 0.875, 0.875]),
    ('min', True, 'no_na', False, [1., 1., 6., 5., 1., 4., 6., 6.]),
    ('min', True, 'no_na', True,
        [0.125, 0.125, 0.75, 0.625, 0.125, 0.5, 0.75, 0.75]),
    ('min', False, 'no_na', False, [3., 3., 6., 1., 3., 2., 6., 6.]),
    ('min', False, 'no_na', True,
        [0.375, 0.375, 0.75, 0.125, 0.375, 0.25, 0.75, 0.75]),
    ('max', True, 'no_na', False, [3., 3., 8., 5., 3., 4., 8., 8.]),
    ('max', True, 'no_na', True,
        [0.375, 0.375, 1., 0.625, 0.375, 0.5, 1., 1.]),
    ('max', False, 'no_na', False, [5., 5., 8., 1., 5., 2., 8., 8.]),
    ('max', False, 'no_na', True,
        [0.625, 0.625, 1., 0.125, 0.625, 0.25, 1., 1.]),
    ('first', True, 'no_na', False, [1., 2., 6., 5., 3., 4., 7., 8.]),
    ('first', True, 'no_na', True,
        [0.125, 0.25, 0.75, 0.625, 0.375, 0.5, 0.875, 1.]),
    ('first', False, 'no_na', False, [3., 4., 6., 1., 5., 2., 7., 8.]),
    ('first', False, 'no_na', True,
        [0.375, 0.5, 0.75, 0.125, 0.625, 0.25, 0.875, 1.]),
    ('dense', True, 'no_na', False, [1., 1., 4., 3., 1., 2., 4., 4.]),
    ('dense', True, 'no_na', True,
        [0.125, 0.125, 0.5, 0.375, 0.125, 0.25, 0.5, 0.5]),
    ('dense', False, 'no_na', False, [3., 3., 4., 1., 3., 2., 4., 4.]),
    ('dense', False, 'no_na', True,
        [0.375, 0.375, 0.5, 0.125, 0.375, 0.25, 0.5, 0.5])
])
def test_rank_args_missing(grps, vals, ties_method, ascending,
                           na_option, pct, exp):
    key = np.repeat(grps, len(vals))
    vals = vals * len(grps)
    df = DataFrame({'key': key, 'val': vals})
    result = df.groupby('key').rank(method=ties_method,
                                    ascending=ascending,
                                    na_option=na_option, pct=pct)

    exp_df = DataFrame(exp * len(grps), columns=['val'])
    tm.assert_frame_equal(result, exp_df)


@pytest.mark.parametrize("pct,exp", [
    (False, [3., 3., 3., 3., 3.]),
    (True, [.6, .6, .6, .6, .6])])
def test_rank_resets_each_group(pct, exp):
    df = DataFrame(
        {'key': ['a', 'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b'],
         'val': [1] * 10}
    )
    result = df.groupby('key').rank(pct=pct)
    exp_df = DataFrame(exp * 2, columns=['val'])
    tm.assert_frame_equal(result, exp_df)


def test_rank_avg_even_vals():
    df = DataFrame({'key': ['a'] * 4, 'val': [1] * 4})
    result = df.groupby('key').rank()
    exp_df = DataFrame([2.5, 2.5, 2.5, 2.5], columns=['val'])
    tm.assert_frame_equal(result, exp_df)


@pytest.mark.parametrize("ties_method", [
    'average', 'min', 'max', 'first', 'dense'])
@pytest.mark.parametrize("ascending", [True, False])
@pytest.mark.parametrize("na_option", ["keep", "top", "bottom"])
@pytest.mark.parametrize("pct", [True, False])
@pytest.mark.parametrize("vals", [
    ['bar', 'bar', 'foo', 'bar', 'baz'],
    ['bar', np.nan, 'foo', np.nan, 'baz']
])
def test_rank_object_raises(ties_method, ascending, na_option,
                            pct, vals):
    df = DataFrame({'key': ['foo'] * 5, 'val': vals})
    with tm.assert_raises_regex(TypeError, "not callable"):
        df.groupby('key').rank(method=ties_method,
                               ascending=ascending,
                               na_option=na_option, pct=pct)
