
import numpy as np
from pandas import date_range, Index
import pandas.util.testing as tm
from pandas.core.reshape.util import cartesian_product


class TestCartesianProduct(object):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result1, result2 = cartesian_product([x, y])
        expected1 = np.array(['A', 'A', 'B', 'B', 'C', 'C'])
        expected2 = np.array([1, 22, 1, 22, 1, 22])
        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    def test_datetimeindex(self):
        # regression test for GitHub issue #6439
        # make sure that the ordering on datetimeindex is consistent
        x = date_range('2000-01-01', periods=2)
        result1, result2 = [Index(y).day for y in cartesian_product([x, x])]
        expected1 = Index([1, 1, 2, 2])
        expected2 = Index([1, 2, 1, 2])
        tm.assert_index_equal(result1, expected1)
        tm.assert_index_equal(result2, expected2)

    def test_empty(self):
        # product of empty factors
        X = [[], [0, 1], []]
        Y = [[], [], ['a', 'b', 'c']]
        for x, y in zip(X, Y):
            expected1 = np.array([], dtype=np.asarray(x).dtype)
            expected2 = np.array([], dtype=np.asarray(y).dtype)
            result1, result2 = cartesian_product([x, y])
            tm.assert_numpy_array_equal(result1, expected1)
            tm.assert_numpy_array_equal(result2, expected2)

        # empty product (empty input):
        result = cartesian_product([])
        expected = []
        assert result == expected

    def test_invalid_input(self):
        invalid_inputs = [1, [1], [1, 2], [[1], 2],
                          'a', ['a'], ['a', 'b'], [['a'], 'b']]
        msg = "Input must be a list-like of list-likes"
        for X in invalid_inputs:
            tm.assert_raises_regex(TypeError, msg, cartesian_product, X=X)
