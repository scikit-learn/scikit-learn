# -*- coding: utf-8 -*-
from __future__ import print_function
from datetime import datetime

import pytest

import numpy as np
import pandas as pd
from pandas import (Index, MultiIndex, CategoricalIndex,
                    DataFrame, Categorical, Series, qcut)
from pandas.util.testing import assert_frame_equal, assert_series_equal
import pandas.util.testing as tm


def cartesian_product_for_groupers(result, args, names):
    """ Reindex to a cartesian production for the groupers,
    preserving the nature (Categorical) of each grouper """

    def f(a):
        if isinstance(a, (CategoricalIndex, Categorical)):
            categories = a.categories
            a = Categorical.from_codes(np.arange(len(categories)),
                                       categories=categories,
                                       ordered=a.ordered)
        return a

    index = pd.MultiIndex.from_product(map(f, args), names=names)
    return result.reindex(index).sort_index()


def test_apply_use_categorical_name(df):
    cats = qcut(df.C, 4)

    def get_stats(group):
        return {'min': group.min(),
                'max': group.max(),
                'count': group.count(),
                'mean': group.mean()}

    result = df.groupby(cats, observed=False).D.apply(get_stats)
    assert result.index.names[0] == 'C'


def test_basic():

    cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                       categories=["a", "b", "c", "d"], ordered=True)
    data = DataFrame({"a": [1, 1, 1, 2, 2, 2, 3, 4, 5], "b": cats})

    exp_index = CategoricalIndex(list('abcd'), name='b', ordered=True)
    expected = DataFrame({'a': [1, 2, 4, np.nan]}, index=exp_index)
    result = data.groupby("b", observed=False).mean()
    tm.assert_frame_equal(result, expected)

    cat1 = Categorical(["a", "a", "b", "b"],
                       categories=["a", "b", "z"], ordered=True)
    cat2 = Categorical(["c", "d", "c", "d"],
                       categories=["c", "d", "y"], ordered=True)
    df = DataFrame({"A": cat1, "B": cat2, "values": [1, 2, 3, 4]})

    # single grouper
    gb = df.groupby("A", observed=False)
    exp_idx = CategoricalIndex(['a', 'b', 'z'], name='A', ordered=True)
    expected = DataFrame({'values': Series([3, 7, 0], index=exp_idx)})
    result = gb.sum()
    tm.assert_frame_equal(result, expected)

    # GH 8623
    x = DataFrame([[1, 'John P. Doe'], [2, 'Jane Dove'],
                   [1, 'John P. Doe']],
                  columns=['person_id', 'person_name'])
    x['person_name'] = Categorical(x.person_name)

    g = x.groupby(['person_id'], observed=False)
    result = g.transform(lambda x: x)
    tm.assert_frame_equal(result, x[['person_name']])

    result = x.drop_duplicates('person_name')
    expected = x.iloc[[0, 1]]
    tm.assert_frame_equal(result, expected)

    def f(x):
        return x.drop_duplicates('person_name').iloc[0]

    result = g.apply(f)
    expected = x.iloc[[0, 1]].copy()
    expected.index = Index([1, 2], name='person_id')
    expected['person_name'] = expected['person_name'].astype('object')
    tm.assert_frame_equal(result, expected)

    # GH 9921
    # Monotonic
    df = DataFrame({"a": [5, 15, 25]})
    c = pd.cut(df.a, bins=[0, 10, 20, 30, 40])

    result = df.a.groupby(c, observed=False).transform(sum)
    tm.assert_series_equal(result, df['a'])

    tm.assert_series_equal(
        df.a.groupby(c, observed=False).transform(lambda xs: np.sum(xs)),
        df['a'])
    tm.assert_frame_equal(
        df.groupby(c, observed=False).transform(sum),
        df[['a']])
    tm.assert_frame_equal(
        df.groupby(c, observed=False).transform(lambda xs: np.max(xs)),
        df[['a']])

    # Filter
    tm.assert_series_equal(
        df.a.groupby(c, observed=False).filter(np.all),
        df['a'])
    tm.assert_frame_equal(
        df.groupby(c, observed=False).filter(np.all),
        df)

    # Non-monotonic
    df = DataFrame({"a": [5, 15, 25, -5]})
    c = pd.cut(df.a, bins=[-10, 0, 10, 20, 30, 40])

    result = df.a.groupby(c, observed=False).transform(sum)
    tm.assert_series_equal(result, df['a'])

    tm.assert_series_equal(
        df.a.groupby(c, observed=False).transform(lambda xs: np.sum(xs)),
        df['a'])
    tm.assert_frame_equal(
        df.groupby(c, observed=False).transform(sum),
        df[['a']])
    tm.assert_frame_equal(
        df.groupby(c, observed=False).transform(lambda xs: np.sum(xs)),
        df[['a']])

    # GH 9603
    df = DataFrame({'a': [1, 0, 0, 0]})
    c = pd.cut(df.a, [0, 1, 2, 3, 4], labels=Categorical(list('abcd')))
    result = df.groupby(c, observed=False).apply(len)

    exp_index = CategoricalIndex(
        c.values.categories, ordered=c.values.ordered)
    expected = Series([1, 0, 0, 0], index=exp_index)
    expected.index.name = 'a'
    tm.assert_series_equal(result, expected)

    # more basic
    levels = ['foo', 'bar', 'baz', 'qux']
    codes = np.random.randint(0, 4, size=100)

    cats = Categorical.from_codes(codes, levels, ordered=True)

    data = DataFrame(np.random.randn(100, 4))

    result = data.groupby(cats, observed=False).mean()

    expected = data.groupby(np.asarray(cats), observed=False).mean()
    exp_idx = CategoricalIndex(levels, categories=cats.categories,
                               ordered=True)
    expected = expected.reindex(exp_idx)

    assert_frame_equal(result, expected)

    grouped = data.groupby(cats, observed=False)
    desc_result = grouped.describe()

    idx = cats.codes.argsort()
    ord_labels = np.asarray(cats).take(idx)
    ord_data = data.take(idx)

    exp_cats = Categorical(ord_labels, ordered=True,
                           categories=['foo', 'bar', 'baz', 'qux'])
    expected = ord_data.groupby(
        exp_cats, sort=False, observed=False).describe()
    assert_frame_equal(desc_result, expected)

    # GH 10460
    expc = Categorical.from_codes(np.arange(4).repeat(8),
                                  levels, ordered=True)
    exp = CategoricalIndex(expc)
    tm.assert_index_equal((desc_result.stack().index
                           .get_level_values(0)), exp)
    exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                 '75%', 'max'] * 4)
    tm.assert_index_equal((desc_result.stack().index
                           .get_level_values(1)), exp)


def test_level_get_group(observed):
    # GH15155
    df = DataFrame(data=np.arange(2, 22, 2),
                   index=MultiIndex(
                       levels=[pd.CategoricalIndex(["a", "b"]), range(10)],
                       labels=[[0] * 5 + [1] * 5, range(10)],
                       names=["Index1", "Index2"]))
    g = df.groupby(level=["Index1"], observed=observed)

    # expected should equal test.loc[["a"]]
    # GH15166
    expected = DataFrame(data=np.arange(2, 12, 2),
                         index=pd.MultiIndex(levels=[pd.CategoricalIndex(
                             ["a", "b"]), range(5)],
        labels=[[0] * 5, range(5)],
        names=["Index1", "Index2"]))
    result = g.get_group('a')

    assert_frame_equal(result, expected)


@pytest.mark.parametrize('ordered', [True, False])
def test_apply(ordered):
    # GH 10138

    dense = Categorical(list('abc'), ordered=ordered)

    # 'b' is in the categories but not in the list
    missing = Categorical(
        list('aaa'), categories=['a', 'b'], ordered=ordered)
    values = np.arange(len(dense))
    df = DataFrame({'missing': missing,
                    'dense': dense,
                    'values': values})
    grouped = df.groupby(['missing', 'dense'], observed=True)

    # missing category 'b' should still exist in the output index
    idx = MultiIndex.from_arrays(
        [missing, dense], names=['missing', 'dense'])
    expected = DataFrame([0, 1, 2.],
                         index=idx,
                         columns=['values'])

    result = grouped.apply(lambda x: np.mean(x))
    assert_frame_equal(result, expected)

    # we coerce back to ints
    expected = expected.astype('int')
    result = grouped.mean()
    assert_frame_equal(result, expected)

    result = grouped.agg(np.mean)
    assert_frame_equal(result, expected)

    # but for transform we should still get back the original index
    idx = MultiIndex.from_arrays([missing, dense],
                                 names=['missing', 'dense'])
    expected = Series(1, index=idx)
    result = grouped.apply(lambda x: 1)
    assert_series_equal(result, expected)


def test_observed(observed):
    # multiple groupers, don't re-expand the output space
    # of the grouper
    # gh-14942 (implement)
    # gh-10132 (back-compat)
    # gh-8138 (back-compat)
    # gh-8869

    cat1 = Categorical(["a", "a", "b", "b"],
                       categories=["a", "b", "z"], ordered=True)
    cat2 = Categorical(["c", "d", "c", "d"],
                       categories=["c", "d", "y"], ordered=True)
    df = DataFrame({"A": cat1, "B": cat2, "values": [1, 2, 3, 4]})
    df['C'] = ['foo', 'bar'] * 2

    # multiple groupers with a non-cat
    gb = df.groupby(['A', 'B', 'C'], observed=observed)
    exp_index = pd.MultiIndex.from_arrays(
        [cat1, cat2, ['foo', 'bar'] * 2],
        names=['A', 'B', 'C'])
    expected = DataFrame({'values': Series(
        [1, 2, 3, 4], index=exp_index)}).sort_index()
    result = gb.sum()
    if not observed:
        expected = cartesian_product_for_groupers(
            expected,
            [cat1, cat2, ['foo', 'bar']],
            list('ABC'))

    tm.assert_frame_equal(result, expected)

    gb = df.groupby(['A', 'B'], observed=observed)
    exp_index = pd.MultiIndex.from_arrays(
        [cat1, cat2],
        names=['A', 'B'])
    expected = DataFrame({'values': [1, 2, 3, 4]},
                         index=exp_index)
    result = gb.sum()
    if not observed:
        expected = cartesian_product_for_groupers(
            expected,
            [cat1, cat2],
            list('AB'))

    tm.assert_frame_equal(result, expected)

    # https://github.com/pandas-dev/pandas/issues/8138
    d = {'cat':
         pd.Categorical(["a", "b", "a", "b"], categories=["a", "b", "c"],
                        ordered=True),
         'ints': [1, 1, 2, 2],
         'val': [10, 20, 30, 40]}
    df = pd.DataFrame(d)

    # Grouping on a single column
    groups_single_key = df.groupby("cat", observed=observed)
    result = groups_single_key.mean()

    exp_index = pd.CategoricalIndex(list('ab'), name="cat",
                                    categories=list('abc'),
                                    ordered=True)
    expected = DataFrame({"ints": [1.5, 1.5], "val": [20., 30]},
                         index=exp_index)
    if not observed:
        index = pd.CategoricalIndex(list('abc'), name="cat",
                                    categories=list('abc'),
                                    ordered=True)
        expected = expected.reindex(index)

    tm.assert_frame_equal(result, expected)

    # Grouping on two columns
    groups_double_key = df.groupby(["cat", "ints"], observed=observed)
    result = groups_double_key.agg('mean')
    expected = DataFrame(
        {"val": [10, 30, 20, 40],
         "cat": pd.Categorical(['a', 'a', 'b', 'b'],
                               categories=['a', 'b', 'c'],
                               ordered=True),
         "ints": [1, 2, 1, 2]}).set_index(["cat", "ints"])
    if not observed:
        expected = cartesian_product_for_groupers(
            expected,
            [df.cat.values, [1, 2]],
            ['cat', 'ints'])

    tm.assert_frame_equal(result, expected)

    # GH 10132
    for key in [('a', 1), ('b', 2), ('b', 1), ('a', 2)]:
        c, i = key
        result = groups_double_key.get_group(key)
        expected = df[(df.cat == c) & (df.ints == i)]
        assert_frame_equal(result, expected)

    # gh-8869
    # with as_index
    d = {'foo': [10, 8, 4, 8, 4, 1, 1], 'bar': [10, 20, 30, 40, 50, 60, 70],
         'baz': ['d', 'c', 'e', 'a', 'a', 'd', 'c']}
    df = pd.DataFrame(d)
    cat = pd.cut(df['foo'], np.linspace(0, 10, 3))
    df['range'] = cat
    groups = df.groupby(['range', 'baz'], as_index=False, observed=observed)
    result = groups.agg('mean')

    groups2 = df.groupby(['range', 'baz'], as_index=True, observed=observed)
    expected = groups2.agg('mean').reset_index()
    tm.assert_frame_equal(result, expected)


def test_observed_codes_remap(observed):
    d = {'C1': [3, 3, 4, 5], 'C2': [1, 2, 3, 4], 'C3': [10, 100, 200, 34]}
    df = pd.DataFrame(d)
    values = pd.cut(df['C1'], [1, 2, 3, 6])
    values.name = "cat"
    groups_double_key = df.groupby([values, 'C2'], observed=observed)

    idx = MultiIndex.from_arrays([values, [1, 2, 3, 4]],
                                 names=["cat", "C2"])
    expected = DataFrame({"C1": [3, 3, 4, 5],
                          "C3": [10, 100, 200, 34]}, index=idx)
    if not observed:
        expected = cartesian_product_for_groupers(
            expected,
            [values.values, [1, 2, 3, 4]],
            ['cat', 'C2'])

    result = groups_double_key.agg('mean')
    tm.assert_frame_equal(result, expected)


def test_observed_perf():
    # we create a cartesian product, so this is
    # non-performant if we don't use observed values
    # gh-14942
    df = DataFrame({
        'cat': np.random.randint(0, 255, size=30000),
        'int_id': np.random.randint(0, 255, size=30000),
        'other_id': np.random.randint(0, 10000, size=30000),
        'foo': 0})
    df['cat'] = df.cat.astype(str).astype('category')

    grouped = df.groupby(['cat', 'int_id', 'other_id'], observed=True)
    result = grouped.count()
    assert result.index.levels[0].nunique() == df.cat.nunique()
    assert result.index.levels[1].nunique() == df.int_id.nunique()
    assert result.index.levels[2].nunique() == df.other_id.nunique()


def test_observed_groups(observed):
    # gh-20583
    # test that we have the appropriate groups

    cat = pd.Categorical(['a', 'c', 'a'], categories=['a', 'b', 'c'])
    df = pd.DataFrame({'cat': cat, 'vals': [1, 2, 3]})
    g = df.groupby('cat', observed=observed)

    result = g.groups
    if observed:
        expected = {'a': Index([0, 2], dtype='int64'),
                    'c': Index([1], dtype='int64')}
    else:
        expected = {'a': Index([0, 2], dtype='int64'),
                    'b': Index([], dtype='int64'),
                    'c': Index([1], dtype='int64')}

    tm.assert_dict_equal(result, expected)


def test_datetime():
    # GH9049: ensure backward compatibility
    levels = pd.date_range('2014-01-01', periods=4)
    codes = np.random.randint(0, 4, size=100)

    cats = Categorical.from_codes(codes, levels, ordered=True)

    data = DataFrame(np.random.randn(100, 4))
    result = data.groupby(cats, observed=False).mean()

    expected = data.groupby(np.asarray(cats), observed=False).mean()
    expected = expected.reindex(levels)
    expected.index = CategoricalIndex(expected.index,
                                      categories=expected.index,
                                      ordered=True)

    assert_frame_equal(result, expected)

    grouped = data.groupby(cats, observed=False)
    desc_result = grouped.describe()

    idx = cats.codes.argsort()
    ord_labels = cats.take_nd(idx)
    ord_data = data.take(idx)
    expected = ord_data.groupby(ord_labels, observed=False).describe()
    assert_frame_equal(desc_result, expected)
    tm.assert_index_equal(desc_result.index, expected.index)
    tm.assert_index_equal(
        desc_result.index.get_level_values(0),
        expected.index.get_level_values(0))

    # GH 10460
    expc = Categorical.from_codes(
        np.arange(4).repeat(8), levels, ordered=True)
    exp = CategoricalIndex(expc)
    tm.assert_index_equal((desc_result.stack().index
                           .get_level_values(0)), exp)
    exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                 '75%', 'max'] * 4)
    tm.assert_index_equal((desc_result.stack().index
                           .get_level_values(1)), exp)


def test_categorical_index():

    s = np.random.RandomState(12345)
    levels = ['foo', 'bar', 'baz', 'qux']
    codes = s.randint(0, 4, size=20)
    cats = Categorical.from_codes(codes, levels, ordered=True)
    df = DataFrame(
        np.repeat(
            np.arange(20), 4).reshape(-1, 4), columns=list('abcd'))
    df['cats'] = cats

    # with a cat index
    result = df.set_index('cats').groupby(level=0, observed=False).sum()
    expected = df[list('abcd')].groupby(cats.codes, observed=False).sum()
    expected.index = CategoricalIndex(
        Categorical.from_codes(
            [0, 1, 2, 3], levels, ordered=True), name='cats')
    assert_frame_equal(result, expected)

    # with a cat column, should produce a cat index
    result = df.groupby('cats', observed=False).sum()
    expected = df[list('abcd')].groupby(cats.codes, observed=False).sum()
    expected.index = CategoricalIndex(
        Categorical.from_codes(
            [0, 1, 2, 3], levels, ordered=True), name='cats')
    assert_frame_equal(result, expected)


def test_describe_categorical_columns():
    # GH 11558
    cats = pd.CategoricalIndex(['qux', 'foo', 'baz', 'bar'],
                               categories=['foo', 'bar', 'baz', 'qux'],
                               ordered=True)
    df = DataFrame(np.random.randn(20, 4), columns=cats)
    result = df.groupby([1, 2, 3, 4] * 5).describe()

    tm.assert_index_equal(result.stack().columns, cats)
    tm.assert_categorical_equal(result.stack().columns.values, cats.values)


def test_unstack_categorical():
    # GH11558 (example is taken from the original issue)
    df = pd.DataFrame({'a': range(10),
                       'medium': ['A', 'B'] * 5,
                       'artist': list('XYXXY') * 2})
    df['medium'] = df['medium'].astype('category')

    gcat = df.groupby(
        ['artist', 'medium'], observed=False)['a'].count().unstack()
    result = gcat.describe()

    exp_columns = pd.CategoricalIndex(['A', 'B'], ordered=False,
                                      name='medium')
    tm.assert_index_equal(result.columns, exp_columns)
    tm.assert_categorical_equal(result.columns.values, exp_columns.values)

    result = gcat['A'] + gcat['B']
    expected = pd.Series([6, 4], index=pd.Index(['X', 'Y'], name='artist'))
    tm.assert_series_equal(result, expected)


def test_bins_unequal_len():
    # GH3011
    series = Series([np.nan, np.nan, 1, 1, 2, 2, 3, 3, 4, 4])
    bins = pd.cut(series.dropna().values, 4)

    # len(bins) != len(series) here
    def f():
        series.groupby(bins).mean()
    pytest.raises(ValueError, f)


def test_as_index():
    # GH13204
    df = DataFrame({'cat': Categorical([1, 2, 2], [1, 2, 3]),
                    'A': [10, 11, 11],
                    'B': [101, 102, 103]})
    result = df.groupby(['cat', 'A'], as_index=False, observed=True).sum()
    expected = DataFrame(
        {'cat': Categorical([1, 2], categories=df.cat.cat.categories),
         'A': [10, 11],
         'B': [101, 205]},
        columns=['cat', 'A', 'B'])
    tm.assert_frame_equal(result, expected)

    # function grouper
    f = lambda r: df.loc[r, 'A']
    result = df.groupby(['cat', f], as_index=False, observed=True).sum()
    expected = DataFrame(
        {'cat': Categorical([1, 2], categories=df.cat.cat.categories),
         'A': [10, 22],
         'B': [101, 205]},
        columns=['cat', 'A', 'B'])
    tm.assert_frame_equal(result, expected)

    # another not in-axis grouper
    s = Series(['a', 'b', 'b'], name='cat2')
    result = df.groupby(['cat', s], as_index=False, observed=True).sum()
    tm.assert_frame_equal(result, expected)

    # GH18872: conflicting names in desired index
    with pytest.raises(ValueError):
        df.groupby(['cat', s.rename('cat')], observed=True).sum()

    # is original index dropped?
    group_columns = ['cat', 'A']
    expected = DataFrame(
        {'cat': Categorical([1, 2], categories=df.cat.cat.categories),
         'A': [10, 11],
         'B': [101, 205]},
        columns=['cat', 'A', 'B'])

    for name in [None, 'X', 'B', 'cat']:
        df.index = Index(list("abc"), name=name)

        if name in group_columns and name in df.index.names:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                result = df.groupby(
                    group_columns, as_index=False, observed=True).sum()

        else:
            result = df.groupby(
                group_columns, as_index=False, observed=True).sum()

        tm.assert_frame_equal(result, expected)


def test_preserve_categories():
    # GH-13179
    categories = list('abc')

    # ordered=True
    df = DataFrame({'A': pd.Categorical(list('ba'),
                                        categories=categories,
                                        ordered=True)})
    index = pd.CategoricalIndex(categories, categories, ordered=True)
    tm.assert_index_equal(
        df.groupby('A', sort=True, observed=False).first().index, index)
    tm.assert_index_equal(
        df.groupby('A', sort=False, observed=False).first().index, index)

    # ordered=False
    df = DataFrame({'A': pd.Categorical(list('ba'),
                                        categories=categories,
                                        ordered=False)})
    sort_index = pd.CategoricalIndex(categories, categories, ordered=False)
    nosort_index = pd.CategoricalIndex(list('bac'), list('bac'),
                                       ordered=False)
    tm.assert_index_equal(
        df.groupby('A', sort=True, observed=False).first().index,
        sort_index)
    tm.assert_index_equal(
        df.groupby('A', sort=False, observed=False).first().index,
        nosort_index)


def test_preserve_categorical_dtype():
    # GH13743, GH13854
    df = DataFrame({'A': [1, 2, 1, 1, 2],
                    'B': [10, 16, 22, 28, 34],
                    'C1': Categorical(list("abaab"),
                                      categories=list("bac"),
                                      ordered=False),
                    'C2': Categorical(list("abaab"),
                                      categories=list("bac"),
                                      ordered=True)})
    # single grouper
    exp_full = DataFrame({'A': [2.0, 1.0, np.nan],
                          'B': [25.0, 20.0, np.nan],
                          'C1': Categorical(list("bac"),
                                            categories=list("bac"),
                                            ordered=False),
                          'C2': Categorical(list("bac"),
                                            categories=list("bac"),
                                            ordered=True)})
    for col in ['C1', 'C2']:
        result1 = df.groupby(by=col, as_index=False, observed=False).mean()
        result2 = df.groupby(
            by=col, as_index=True, observed=False).mean().reset_index()
        expected = exp_full.reindex(columns=result1.columns)
        tm.assert_frame_equal(result1, expected)
        tm.assert_frame_equal(result2, expected)


def test_categorical_no_compress():
    data = Series(np.random.randn(9))

    codes = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    cats = Categorical.from_codes(codes, [0, 1, 2], ordered=True)

    result = data.groupby(cats, observed=False).mean()
    exp = data.groupby(codes, observed=False).mean()

    exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                 ordered=cats.ordered)
    assert_series_equal(result, exp)

    codes = np.array([0, 0, 0, 1, 1, 1, 3, 3, 3])
    cats = Categorical.from_codes(codes, [0, 1, 2, 3], ordered=True)

    result = data.groupby(cats, observed=False).mean()
    exp = data.groupby(codes, observed=False).mean().reindex(cats.categories)
    exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                 ordered=cats.ordered)
    assert_series_equal(result, exp)

    cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                       categories=["a", "b", "c", "d"], ordered=True)
    data = DataFrame({"a": [1, 1, 1, 2, 2, 2, 3, 4, 5], "b": cats})

    result = data.groupby("b", observed=False).mean()
    result = result["a"].values
    exp = np.array([1, 2, 4, np.nan])
    tm.assert_numpy_array_equal(result, exp)


def test_sort():

    # http://stackoverflow.com/questions/23814368/sorting-pandas-categorical-labels-after-groupby  # noqa: flake8
    # This should result in a properly sorted Series so that the plot
    # has a sorted x axis
    # self.cat.groupby(['value_group'])['value_group'].count().plot(kind='bar')

    df = DataFrame({'value': np.random.randint(0, 10000, 100)})
    labels = ["{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500)]
    cat_labels = Categorical(labels, labels)

    df = df.sort_values(by=['value'], ascending=True)
    df['value_group'] = pd.cut(df.value, range(0, 10500, 500),
                               right=False, labels=cat_labels)

    res = df.groupby(['value_group'], observed=False)['value_group'].count()
    exp = res[sorted(res.index, key=lambda x: float(x.split()[0]))]
    exp.index = CategoricalIndex(exp.index, name=exp.index.name)
    tm.assert_series_equal(res, exp)


def test_sort2():
    # dataframe groupby sort was being ignored # GH 8868
    df = DataFrame([['(7.5, 10]', 10, 10],
                    ['(7.5, 10]', 8, 20],
                    ['(2.5, 5]', 5, 30],
                    ['(5, 7.5]', 6, 40],
                    ['(2.5, 5]', 4, 50],
                    ['(0, 2.5]', 1, 60],
                    ['(5, 7.5]', 7, 70]], columns=['range', 'foo', 'bar'])
    df['range'] = Categorical(df['range'], ordered=True)
    index = CategoricalIndex(['(0, 2.5]', '(2.5, 5]', '(5, 7.5]',
                              '(7.5, 10]'], name='range', ordered=True)
    expected_sort = DataFrame([[1, 60], [5, 30], [6, 40], [10, 10]],
                              columns=['foo', 'bar'], index=index)

    col = 'range'
    result_sort = df.groupby(col, sort=True, observed=False).first()
    assert_frame_equal(result_sort, expected_sort)

    # when categories is ordered, group is ordered by category's order
    expected_sort = result_sort
    result_sort = df.groupby(col, sort=False, observed=False).first()
    assert_frame_equal(result_sort, expected_sort)

    df['range'] = Categorical(df['range'], ordered=False)
    index = CategoricalIndex(['(0, 2.5]', '(2.5, 5]', '(5, 7.5]',
                              '(7.5, 10]'], name='range')
    expected_sort = DataFrame([[1, 60], [5, 30], [6, 40], [10, 10]],
                              columns=['foo', 'bar'], index=index)

    index = CategoricalIndex(['(7.5, 10]', '(2.5, 5]', '(5, 7.5]',
                              '(0, 2.5]'],
                             categories=['(7.5, 10]', '(2.5, 5]',
                                         '(5, 7.5]', '(0, 2.5]'],
                             name='range')
    expected_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                                index=index, columns=['foo', 'bar'])

    col = 'range'

    # this is an unordered categorical, but we allow this ####
    result_sort = df.groupby(col, sort=True, observed=False).first()
    assert_frame_equal(result_sort, expected_sort)

    result_nosort = df.groupby(col, sort=False, observed=False).first()
    assert_frame_equal(result_nosort, expected_nosort)


def test_sort_datetimelike():
    # GH10505

    # use same data as test_groupby_sort_categorical, which category is
    # corresponding to datetime.month
    df = DataFrame({'dt': [datetime(2011, 7, 1), datetime(2011, 7, 1),
                           datetime(2011, 2, 1), datetime(2011, 5, 1),
                           datetime(2011, 2, 1), datetime(2011, 1, 1),
                           datetime(2011, 5, 1)],
                    'foo': [10, 8, 5, 6, 4, 1, 7],
                    'bar': [10, 20, 30, 40, 50, 60, 70]},
                   columns=['dt', 'foo', 'bar'])

    # ordered=True
    df['dt'] = Categorical(df['dt'], ordered=True)
    index = [datetime(2011, 1, 1), datetime(2011, 2, 1),
             datetime(2011, 5, 1), datetime(2011, 7, 1)]
    result_sort = DataFrame(
        [[1, 60], [5, 30], [6, 40], [10, 10]], columns=['foo', 'bar'])
    result_sort.index = CategoricalIndex(index, name='dt', ordered=True)

    index = [datetime(2011, 7, 1), datetime(2011, 2, 1),
             datetime(2011, 5, 1), datetime(2011, 1, 1)]
    result_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                              columns=['foo', 'bar'])
    result_nosort.index = CategoricalIndex(index, categories=index,
                                           name='dt', ordered=True)

    col = 'dt'
    assert_frame_equal(
        result_sort, df.groupby(col, sort=True, observed=False).first())

    # when categories is ordered, group is ordered by category's order
    assert_frame_equal(
        result_sort, df.groupby(col, sort=False, observed=False).first())

    # ordered = False
    df['dt'] = Categorical(df['dt'], ordered=False)
    index = [datetime(2011, 1, 1), datetime(2011, 2, 1),
             datetime(2011, 5, 1), datetime(2011, 7, 1)]
    result_sort = DataFrame(
        [[1, 60], [5, 30], [6, 40], [10, 10]], columns=['foo', 'bar'])
    result_sort.index = CategoricalIndex(index, name='dt')

    index = [datetime(2011, 7, 1), datetime(2011, 2, 1),
             datetime(2011, 5, 1), datetime(2011, 1, 1)]
    result_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                              columns=['foo', 'bar'])
    result_nosort.index = CategoricalIndex(index, categories=index,
                                           name='dt')

    col = 'dt'
    assert_frame_equal(
        result_sort, df.groupby(col, sort=True, observed=False).first())
    assert_frame_equal(
        result_nosort, df.groupby(col, sort=False, observed=False).first())


def test_empty_sum():
    # https://github.com/pandas-dev/pandas/issues/18678
    df = pd.DataFrame({"A": pd.Categorical(['a', 'a', 'b'],
                                           categories=['a', 'b', 'c']),
                       'B': [1, 2, 1]})
    expected_idx = pd.CategoricalIndex(['a', 'b', 'c'], name='A')

    # 0 by default
    result = df.groupby("A", observed=False).B.sum()
    expected = pd.Series([3, 1, 0], expected_idx, name='B')
    tm.assert_series_equal(result, expected)

    # min_count=0
    result = df.groupby("A", observed=False).B.sum(min_count=0)
    expected = pd.Series([3, 1, 0], expected_idx, name='B')
    tm.assert_series_equal(result, expected)

    # min_count=1
    result = df.groupby("A", observed=False).B.sum(min_count=1)
    expected = pd.Series([3, 1, np.nan], expected_idx, name='B')
    tm.assert_series_equal(result, expected)

    # min_count>1
    result = df.groupby("A", observed=False).B.sum(min_count=2)
    expected = pd.Series([3, np.nan, np.nan], expected_idx, name='B')
    tm.assert_series_equal(result, expected)


def test_empty_prod():
    # https://github.com/pandas-dev/pandas/issues/18678
    df = pd.DataFrame({"A": pd.Categorical(['a', 'a', 'b'],
                                           categories=['a', 'b', 'c']),
                       'B': [1, 2, 1]})

    expected_idx = pd.CategoricalIndex(['a', 'b', 'c'], name='A')

    # 1 by default
    result = df.groupby("A", observed=False).B.prod()
    expected = pd.Series([2, 1, 1], expected_idx, name='B')
    tm.assert_series_equal(result, expected)

    # min_count=0
    result = df.groupby("A", observed=False).B.prod(min_count=0)
    expected = pd.Series([2, 1, 1], expected_idx, name='B')
    tm.assert_series_equal(result, expected)

    # min_count=1
    result = df.groupby("A", observed=False).B.prod(min_count=1)
    expected = pd.Series([2, 1, np.nan], expected_idx, name='B')
    tm.assert_series_equal(result, expected)
