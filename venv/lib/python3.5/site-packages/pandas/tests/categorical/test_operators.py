# -*- coding: utf-8 -*-

import pytest

import pandas as pd
import numpy as np

import pandas.util.testing as tm
from pandas import Categorical, Series, DataFrame, date_range
from pandas.tests.categorical.common import TestCategorical


class TestCategoricalOpsWithFactor(TestCategorical):

    def test_categories_none_comparisons(self):
        factor = Categorical(['a', 'b', 'b', 'a',
                              'a', 'c', 'c', 'c'], ordered=True)
        tm.assert_categorical_equal(factor, self.factor)

    def test_comparisons(self):

        result = self.factor[self.factor == 'a']
        expected = self.factor[np.asarray(self.factor) == 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor != 'a']
        expected = self.factor[np.asarray(self.factor) != 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor < 'c']
        expected = self.factor[np.asarray(self.factor) < 'c']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor > 'a']
        expected = self.factor[np.asarray(self.factor) > 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor >= 'b']
        expected = self.factor[np.asarray(self.factor) >= 'b']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor <= 'b']
        expected = self.factor[np.asarray(self.factor) <= 'b']
        tm.assert_categorical_equal(result, expected)

        n = len(self.factor)

        other = self.factor[np.random.permutation(n)]
        result = self.factor == other
        expected = np.asarray(self.factor) == np.asarray(other)
        tm.assert_numpy_array_equal(result, expected)

        result = self.factor == 'd'
        expected = np.repeat(False, len(self.factor))
        tm.assert_numpy_array_equal(result, expected)

        # comparisons with categoricals
        cat_rev = Categorical(
            ["a", "b", "c"], categories=["c", "b", "a"], ordered=True)
        cat_rev_base = Categorical(
            ["b", "b", "b"], categories=["c", "b", "a"], ordered=True)
        cat = Categorical(["a", "b", "c"], ordered=True)
        cat_base = Categorical(
            ["b", "b", "b"], categories=cat.categories, ordered=True)

        # comparisons need to take categories ordering into account
        res_rev = cat_rev > cat_rev_base
        exp_rev = np.array([True, False, False])
        tm.assert_numpy_array_equal(res_rev, exp_rev)

        res_rev = cat_rev < cat_rev_base
        exp_rev = np.array([False, False, True])
        tm.assert_numpy_array_equal(res_rev, exp_rev)

        res = cat > cat_base
        exp = np.array([False, False, True])
        tm.assert_numpy_array_equal(res, exp)

        # Only categories with same categories can be compared
        def f():
            cat > cat_rev

        pytest.raises(TypeError, f)

        cat_rev_base2 = Categorical(
            ["b", "b", "b"], categories=["c", "b", "a", "d"])

        def f():
            cat_rev > cat_rev_base2

        pytest.raises(TypeError, f)

        # Only categories with same ordering information can be compared
        cat_unorderd = cat.set_ordered(False)
        assert not (cat > cat).any()

        def f():
            cat > cat_unorderd

        pytest.raises(TypeError, f)

        # comparison (in both directions) with Series will raise
        s = Series(["b", "b", "b"])
        pytest.raises(TypeError, lambda: cat > s)
        pytest.raises(TypeError, lambda: cat_rev > s)
        pytest.raises(TypeError, lambda: s < cat)
        pytest.raises(TypeError, lambda: s < cat_rev)

        # comparison with numpy.array will raise in both direction, but only on
        # newer numpy versions
        a = np.array(["b", "b", "b"])
        pytest.raises(TypeError, lambda: cat > a)
        pytest.raises(TypeError, lambda: cat_rev > a)

        # Make sure that unequal comparison take the categories order in
        # account
        cat_rev = Categorical(
            list("abc"), categories=list("cba"), ordered=True)
        exp = np.array([True, False, False])
        res = cat_rev > "b"
        tm.assert_numpy_array_equal(res, exp)


class TestCategoricalOps(object):

    def test_datetime_categorical_comparison(self):
        dt_cat = Categorical(date_range('2014-01-01', periods=3), ordered=True)
        tm.assert_numpy_array_equal(dt_cat > dt_cat[0],
                                    np.array([False, True, True]))
        tm.assert_numpy_array_equal(dt_cat[0] < dt_cat,
                                    np.array([False, True, True]))

    def test_reflected_comparison_with_scalars(self):
        # GH8658
        cat = Categorical([1, 2, 3], ordered=True)
        tm.assert_numpy_array_equal(cat > cat[0],
                                    np.array([False, True, True]))
        tm.assert_numpy_array_equal(cat[0] < cat,
                                    np.array([False, True, True]))

    def test_comparison_with_unknown_scalars(self):
        # https://github.com/pandas-dev/pandas/issues/9836#issuecomment-92123057
        # and following comparisons with scalars not in categories should raise
        # for unequal comps, but not for equal/not equal
        cat = Categorical([1, 2, 3], ordered=True)

        pytest.raises(TypeError, lambda: cat < 4)
        pytest.raises(TypeError, lambda: cat > 4)
        pytest.raises(TypeError, lambda: 4 < cat)
        pytest.raises(TypeError, lambda: 4 > cat)

        tm.assert_numpy_array_equal(cat == 4,
                                    np.array([False, False, False]))
        tm.assert_numpy_array_equal(cat != 4,
                                    np.array([True, True, True]))

    @pytest.mark.parametrize('data,reverse,base', [
        (list("abc"), list("cba"), list("bbb")),
        ([1, 2, 3], [3, 2, 1], [2, 2, 2])]
    )
    def test_comparisons(self, data, reverse, base):
        cat_rev = Series(
            Categorical(data, categories=reverse, ordered=True))
        cat_rev_base = Series(
            Categorical(base, categories=reverse, ordered=True))
        cat = Series(Categorical(data, ordered=True))
        cat_base = Series(
            Categorical(base, categories=cat.cat.categories, ordered=True))
        s = Series(base)
        a = np.array(base)

        # comparisons need to take categories ordering into account
        res_rev = cat_rev > cat_rev_base
        exp_rev = Series([True, False, False])
        tm.assert_series_equal(res_rev, exp_rev)

        res_rev = cat_rev < cat_rev_base
        exp_rev = Series([False, False, True])
        tm.assert_series_equal(res_rev, exp_rev)

        res = cat > cat_base
        exp = Series([False, False, True])
        tm.assert_series_equal(res, exp)

        scalar = base[1]
        res = cat > scalar
        exp = Series([False, False, True])
        exp2 = cat.values > scalar
        tm.assert_series_equal(res, exp)
        tm.assert_numpy_array_equal(res.values, exp2)
        res_rev = cat_rev > scalar
        exp_rev = Series([True, False, False])
        exp_rev2 = cat_rev.values > scalar
        tm.assert_series_equal(res_rev, exp_rev)
        tm.assert_numpy_array_equal(res_rev.values, exp_rev2)

        # Only categories with same categories can be compared
        def f():
            cat > cat_rev

        pytest.raises(TypeError, f)

        # categorical cannot be compared to Series or numpy array, and also
        # not the other way around
        pytest.raises(TypeError, lambda: cat > s)
        pytest.raises(TypeError, lambda: cat_rev > s)
        pytest.raises(TypeError, lambda: cat > a)
        pytest.raises(TypeError, lambda: cat_rev > a)

        pytest.raises(TypeError, lambda: s < cat)
        pytest.raises(TypeError, lambda: s < cat_rev)

        pytest.raises(TypeError, lambda: a < cat)
        pytest.raises(TypeError, lambda: a < cat_rev)

    @pytest.mark.parametrize('ctor', [
        lambda *args, **kwargs: Categorical(*args, **kwargs),
        lambda *args, **kwargs: Series(Categorical(*args, **kwargs)),
    ])
    def test_unordered_different_order_equal(self, ctor):
        # https://github.com/pandas-dev/pandas/issues/16014
        c1 = ctor(['a', 'b'], categories=['a', 'b'], ordered=False)
        c2 = ctor(['a', 'b'], categories=['b', 'a'], ordered=False)
        assert (c1 == c2).all()

        c1 = ctor(['a', 'b'], categories=['a', 'b'], ordered=False)
        c2 = ctor(['b', 'a'], categories=['b', 'a'], ordered=False)
        assert (c1 != c2).all()

        c1 = ctor(['a', 'a'], categories=['a', 'b'], ordered=False)
        c2 = ctor(['b', 'b'], categories=['b', 'a'], ordered=False)
        assert (c1 != c2).all()

        c1 = ctor(['a', 'a'], categories=['a', 'b'], ordered=False)
        c2 = ctor(['a', 'b'], categories=['b', 'a'], ordered=False)
        result = c1 == c2
        tm.assert_numpy_array_equal(np.array(result), np.array([True, False]))

    def test_unordered_different_categories_raises(self):
        c1 = Categorical(['a', 'b'], categories=['a', 'b'], ordered=False)
        c2 = Categorical(['a', 'c'], categories=['c', 'a'], ordered=False)
        with tm.assert_raises_regex(TypeError,
                                    "Categoricals can only be compared"):
            c1 == c2

    def test_compare_different_lengths(self):
        c1 = Categorical([], categories=['a', 'b'])
        c2 = Categorical([], categories=['a'])
        msg = "Categories are different lengths"
        with tm.assert_raises_regex(TypeError, msg):
            c1 == c2

    def test_compare_unordered_different_order(self):
        # https://github.com/pandas-dev/pandas/issues/16603#issuecomment-
        # 349290078
        a = pd.Categorical(['a'], categories=['a', 'b'])
        b = pd.Categorical(['b'], categories=['b', 'a'])
        assert not a.equals(b)

    def test_numeric_like_ops(self):

        df = DataFrame({'value': np.random.randint(0, 10000, 100)})
        labels = ["{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500)]
        cat_labels = Categorical(labels, labels)

        df = df.sort_values(by=['value'], ascending=True)
        df['value_group'] = pd.cut(df.value, range(0, 10500, 500),
                                   right=False, labels=cat_labels)

        # numeric ops should not succeed
        for op in ['__add__', '__sub__', '__mul__', '__truediv__']:
            pytest.raises(TypeError,
                          lambda: getattr(df, op)(df))

        # reduction ops should not succeed (unless specifically defined, e.g.
        # min/max)
        s = df['value_group']
        for op in ['kurt', 'skew', 'var', 'std', 'mean', 'sum', 'median']:
            pytest.raises(TypeError,
                          lambda: getattr(s, op)(numeric_only=False))

        # mad technically works because it takes always the numeric data

        # numpy ops
        s = Series(Categorical([1, 2, 3, 4]))
        pytest.raises(TypeError, lambda: np.sum(s))

        # numeric ops on a Series
        for op in ['__add__', '__sub__', '__mul__', '__truediv__']:
            pytest.raises(TypeError, lambda: getattr(s, op)(2))

        # invalid ufunc
        pytest.raises(TypeError, lambda: np.log(s))
