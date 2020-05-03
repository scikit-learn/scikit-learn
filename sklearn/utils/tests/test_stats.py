import numpy as np
from sklearn.utils.testing import assert_equal, assert_raises
from sklearn.utils.stats import _weighted_percentile


def test_weighted_percentile_negative_weights_raises():
    weight = np.array([1, -1])
    data = np.array([0, 1])
    assert_raises(ValueError, _weighted_percentile, data, weight)


def test_weighted_percentile_negative_percentile_raises():
    weight = np.array([1, -1])
    data = np.array([0, 1])
    percentile = -50
    assert_raises(ValueError, _weighted_percentile, data, weight,
                  percentile=percentile)


def test_weighted_percentile_no_weights():
    weight = np.array([0, 0])
    data = np.array([0, 1])
    percentile = 50
    expected = 0
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_median_interpolated_list():
    weight = [1, 1]
    data = [0, 1]
    percentile = 50
    expected = 0.5
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_median_interpolated_tuple():
    weight = (1, 1)
    data = (0, 1)
    percentile = 50
    expected = 0.5
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_median_interpolated():
    weight = np.array([1, 1])
    data = np.array([0, 1])
    percentile = 50
    expected = 0.5
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_median_regular():
    weight = np.array([1, 1, 1])
    data = np.array([0, 1, 2])
    percentile = 50
    expected = 1.0
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_0_regular():
    weight = np.array([1, 1, 1])
    data = np.array([0, 1, 2])
    percentile = 0
    expected = 0.0
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_90_regular():
    weight = np.array([1, 1, 1])
    data = np.array([0, 1, 2])
    percentile = 90
    expected = 2.0
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_70_interpolated():
    weight = np.array([1, 1, 1, 1])
    data = np.arange(0, 4, 1)
    percentile = 70
    expected = 2.3
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)


def test_weighted_percentile_70_mixed_weights():
    weight = np.array([1, 0, 1, 1])
    data = np.arange(0, 4, 1)
    percentile = 50
    expected = 2.0
    actual = _weighted_percentile(data, weight, percentile=percentile)
    assert_equal(expected, actual)
