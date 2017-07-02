from __future__ import absolute_import

import numpy as np
from six.moves import range
import warnings

from sklearn.preprocessing import KBinsDiscretizer
from sklearn.utils.testing import (
    assert_array_equal,
    assert_raises,
    assert_raise_message,
    assert_warns_message
)

X = np.array([[-2, 1.5, -4, -1],
              [-1, 2.5, -3, -0.5],
              [0, 3.5, -2, 0.5],
              [1, 4.5, -1, 2]])


def test_fit_transform():
    est = KBinsDiscretizer(n_bins=3).fit(X)
    expected = [[0, 0, 0, 0],
                [1, 1, 1, 0],
                [2, 2, 2, 1],
                [2, 2, 2, 2]]
    assert_array_equal(expected, est.transform(X))


def test_invalid_n_bins():
    est = KBinsDiscretizer(n_bins=1)
    assert_raise_message(ValueError, "KBinsDiscretizer received an invalid "
                         "number of bins. Received 1, expected at least 2.",
                         est.fit_transform, X)

    est = KBinsDiscretizer(n_bins=1.1)
    assert_raise_message(ValueError, "KBinsDiscretizer received an invalid "
                         "n_bins type. Received float, expected int.",
                         est.fit_transform, X)


def test_invalid_n_bins_array():
    # Bad shape
    n_bins = np.ones((2, 4)) * 2
    est = KBinsDiscretizer(n_bins=n_bins)
    assert_raise_message(ValueError,
                         "n_bins must be a scalar or array of shape "
                         "(n_features,).", est.fit_transform, X)

    # Incorrect number of features
    n_bins = [1, 2, 2]
    est = KBinsDiscretizer(n_bins=n_bins)
    assert_raise_message(ValueError,
                         "n_bins must be a scalar or array of shape "
                         "(n_features,).", est.fit_transform, X)

    # Bad bin values
    n_bins = [1, 2, 2, 1]
    est = KBinsDiscretizer(n_bins=n_bins)
    assert_raise_message(ValueError,
                         "KBinsDiscretizer received an invalid number of bins "
                         "at indices 0, 3. Number of bins must be at least 2, "
                         "and must be an int.",
                         est.fit_transform, X)

    # Float bin values
    n_bins = [2.1, 2, 2.1, 2]
    est = KBinsDiscretizer(n_bins=n_bins)
    assert_raise_message(ValueError,
                         "KBinsDiscretizer received an invalid number of bins "
                         "at indices 0, 2. Number of bins must be at least 2, "
                         "and must be an int.",
                         est.fit_transform, X)


def test_fit_transform_n_bins_array():
    est = KBinsDiscretizer(n_bins=[2, 3, 3, 3]).fit(X)
    expected = [[0, 0, 0, 0],
                [0, 1, 1, 0],
                [1, 2, 2, 1],
                [1, 2, 2, 2]]
    assert_array_equal(expected, est.transform(X))


def test_invalid_n_features():
    est = KBinsDiscretizer(n_bins=3).fit(X)
    bad_X = np.arange(25).reshape(5, -1)
    assert_raise_message(ValueError,
                         "Incorrect number of features. Expecting 4, "
                         "received 5", est.transform, bad_X)


def test_ignored_transform():
    # Feature at col_idx=1 should not change
    est = KBinsDiscretizer(n_bins=3, ignored_features=[1]).fit(X)

    expected = [[0., 1.5, 0., 0.],
                [1., 2.5, 1., 0.],
                [2., 3.5, 2., 1.],
                [2., 4.5, 2., 2.]]

    assert_array_equal(expected, est.transform(X))


def test_ignored_invalid():
    # Duplicate column
    est = KBinsDiscretizer(ignored_features=[1, 1])
    assert_raise_message(ValueError, "Duplicate ignored column indices found.",
                         est.fit, X)

    invalid_ignored = [
        [-1],  # Invalid index
        [4],  # Invalid index
        ['a'],  # Not an integer index
        [4.5],  # Not an integer index
        [[1, 2], [3, 4]]  # Invalid shape
    ]

    for invalid in invalid_ignored:
        est = KBinsDiscretizer(ignored_features=invalid)
        assert_raises(ValueError, est.fit, X)


def test_same_min_max():
    warnings.simplefilter("always")
    X = np.array([[1, -2],
                  [1, -1],
                  [1, 0],
                  [1, 1]])
    est = assert_warns_message(UserWarning,
                               "Features 0 are constant and will be replaced "
                               "with 0.", KBinsDiscretizer(n_bins=3).fit, X)
    Xt = est.transform(X)

    expected = [[0, 0],
                [0, 1],
                [0, 2],
                [0, 2]]
    assert_array_equal(expected, Xt)


def test_transform_1d_behavior():
    X = np.arange(4)
    est = KBinsDiscretizer(n_bins=2)
    assert_raises(ValueError, est.fit, X)

    est = KBinsDiscretizer(n_bins=2)
    est.fit(X.reshape(-1, 1))
    assert_raises(ValueError, est.transform, X)


def test_inverse_transform_with_ignored():
    est = KBinsDiscretizer(n_bins=[2, 3, 0, 3], ignored_features=[1, 2]).fit(X)
    Xt = [[0, 1, -4.5, 0],
          [0, 2, -3.5, 0],
          [1, 3, -2.5, 1],
          [1, 3, -1.5, 2]]

    Xinv = est.inverse_transform(Xt)
    expected = [[-1.25, 1, -4.5, -0.5],
                [-1.25, 2, -3.5, -0.5],
                [0.25, 3, -2.5, 0.5],
                [0.25, 3, -1.5, 1.5]]

    assert_array_equal(expected, Xinv)


def test_numeric_stability():
    X_init = np.array([2., 4., 6., 8., 10.]).reshape(-1, 1)
    Xt_expected = np.array([0, 0, 1, 1, 1]).reshape(-1, 1)

    # Test up to discretizing nano units
    for i in range(1, 9):
        X = X_init / 10**i
        Xt = KBinsDiscretizer(n_bins=2).fit_transform(X)
        assert_array_equal(Xt_expected, Xt)
