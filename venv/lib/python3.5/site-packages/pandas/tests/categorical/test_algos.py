import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm


@pytest.mark.parametrize('ordered', [True, False])
@pytest.mark.parametrize('categories', [
    ['b', 'a', 'c'],
    ['a', 'b', 'c', 'd'],
])
def test_factorize(categories, ordered):
    cat = pd.Categorical(['b', 'b', 'a', 'c', None],
                         categories=categories,
                         ordered=ordered)
    labels, uniques = pd.factorize(cat)
    expected_labels = np.array([0, 0, 1, 2, -1], dtype=np.intp)
    expected_uniques = pd.Categorical(['b', 'a', 'c'],
                                      categories=categories,
                                      ordered=ordered)

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort():
    cat = pd.Categorical(['b', 'b', None, 'a'])
    labels, uniques = pd.factorize(cat, sort=True)
    expected_labels = np.array([1, 1, -1, 0], dtype=np.intp)
    expected_uniques = pd.Categorical(['a', 'b'])

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort_ordered():
    cat = pd.Categorical(['b', 'b', None, 'a'],
                         categories=['c', 'b', 'a'],
                         ordered=True)

    labels, uniques = pd.factorize(cat, sort=True)
    expected_labels = np.array([0, 0, -1, 1], dtype=np.intp)
    expected_uniques = pd.Categorical(['b', 'a'],
                                      categories=['c', 'b', 'a'],
                                      ordered=True)

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_isin_cats():
    # GH2003
    cat = pd.Categorical(["a", "b", np.nan])

    result = cat.isin(["a", np.nan])
    expected = np.array([True, False, True], dtype=bool)
    tm.assert_numpy_array_equal(expected, result)

    result = cat.isin(["a", "c"])
    expected = np.array([True, False, False], dtype=bool)
    tm.assert_numpy_array_equal(expected, result)


@pytest.mark.parametrize("empty", [[], pd.Series(), np.array([])])
def test_isin_empty(empty):
    s = pd.Categorical(["a", "b"])
    expected = np.array([False, False], dtype=bool)

    result = s.isin(empty)
    tm.assert_numpy_array_equal(expected, result)


class TestTake(object):
    # https://github.com/pandas-dev/pandas/issues/20664

    def test_take_warns(self):
        cat = pd.Categorical(['a', 'b'])
        with tm.assert_produces_warning(FutureWarning):
            cat.take([0, -1])

    def test_take_positive_no_warning(self):
        cat = pd.Categorical(['a', 'b'])
        with tm.assert_produces_warning(None):
            cat.take([0, 0])

    def test_take_bounds(self, allow_fill):
        # https://github.com/pandas-dev/pandas/issues/20664
        cat = pd.Categorical(['a', 'b', 'a'])
        with pytest.raises(IndexError):
            cat.take([4, 5], allow_fill=allow_fill)

    def test_take_empty(self, allow_fill):
        # https://github.com/pandas-dev/pandas/issues/20664
        cat = pd.Categorical([], categories=['a', 'b'])
        with pytest.raises(IndexError):
            cat.take([0], allow_fill=allow_fill)

    def test_positional_take(self, ordered):
        cat = pd.Categorical(['a', 'a', 'b', 'b'], categories=['b', 'a'],
                             ordered=ordered)
        result = cat.take([0, 1, 2], allow_fill=False)
        expected = pd.Categorical(['a', 'a', 'b'], categories=cat.categories,
                                  ordered=ordered)
        tm.assert_categorical_equal(result, expected)

    def test_positional_take_unobserved(self, ordered):
        cat = pd.Categorical(['a', 'b'], categories=['a', 'b', 'c'],
                             ordered=ordered)
        result = cat.take([1, 0], allow_fill=False)
        expected = pd.Categorical(['b', 'a'], categories=cat.categories,
                                  ordered=ordered)
        tm.assert_categorical_equal(result, expected)
