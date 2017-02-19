from __future__ import absolute_import

import numpy as np
import warnings

from sklearn.preprocessing import KBinsDiscretizer
from sklearn.utils.testing import (
    assert_array_equal,
    assert_raises,
    assert_warns
)

X = np.array([[-2, 1, -4, -1],
              [-1, 2, -3, -0.5],
              [0, 3, -2, 0.5],
              [1, 4, -1, 2]])


def test_fit_transform():
    dis = KBinsDiscretizer(n_bins=3).fit(X)
    expected = [[0, 0, 0, 0],
                [1, 1, 1, 0],
                [2, 2, 2, 1],
                [2, 2, 2, 2]]
    assert_array_equal(expected, dis.transform(X))


def test_invalid_n_bins():
    dis = KBinsDiscretizer(n_bins=1)
    assert_raises(ValueError, dis.fit_transform, X)


def test_invalid_n_bins_array():

    # Bad shape
    n_bins = np.ones((2, 4)) * 2
    dis = KBinsDiscretizer(n_bins=n_bins)
    assert_raises(ValueError, dis.fit_transform, X)

    # Bad values
    n_bins = [1, 2, 2, 2]
    dis = KBinsDiscretizer(n_bins=n_bins)
    assert_raises(ValueError, dis.fit_transform, X)


def test_fit_transform_n_bins_array():
    dis = KBinsDiscretizer(n_bins=[2, 3, 3, 3]).fit(X)
    expected = [[0, 0, 0, 0],
                [0, 1, 1, 0],
                [1, 2, 2, 1],
                [1, 2, 2, 2]]
    assert_array_equal(expected, dis.transform(X))


def test_invalid_n_features():
    dis = KBinsDiscretizer(n_bins=3).fit(X)
    bad_X = np.arange(25).reshape(5, -1)
    assert_raises(ValueError, dis.transform, bad_X)


def test_categorical_transform():
    # Feature at col_idx=1 should not change
    dis = KBinsDiscretizer(n_bins=3, categorical_features=[1]).fit(X)

    expected = [[0., 1, 0., 0.],
                [1., 2, 1., 0.],
                [2., 3, 2., 1.],
                [2., 4, 2., 2.]]

    assert_array_equal(expected, dis.transform(X))


def test_categorical_invalid():
    invalid_categorical = [
        [1, 1],  # Duplicate column
        [-1],  # Invalid index
        [4],  # Invalid index
        ['a'],  # Not an integer index
        [4.5],  # Not an integer index
        [[1, 2], [3, 4]]  # Invalid shape
    ]

    for invalid in invalid_categorical:
        dis = KBinsDiscretizer(categorical_features=invalid)
        assert_raises(ValueError, dis.fit, X)


def test_min_max_same():
    warnings.simplefilter("always")
    X = np.array([[1, -2],
                  [1, -1],
                  [1, 0],
                  [1, 1]])
    dis = KBinsDiscretizer(n_bins=3).fit(X)
    X_t = assert_warns(UserWarning, dis.transform, X)

    expected = [[0, 0],
                [0, 1],
                [0, 2],
                [0, 2]]
    assert_array_equal(expected, X_t)


def test_transform_1d_behavior():
    X = np.arange(4)
    dis = KBinsDiscretizer(n_bins=2)
    assert_raises(ValueError, dis.fit, X)

    dis = KBinsDiscretizer(n_bins=2)
    dis.fit(X.reshape(-1, 1))
    assert_raises(ValueError, dis.transform, X)
