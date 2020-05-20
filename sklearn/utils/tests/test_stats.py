import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_almost_equal

from sklearn.utils.stats import _weighted_percentile


def test_weighted_percentile():
    y = np.empty(102, dtype=np.float64)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert score == 1

    # Check for multioutput
    y_multi = np.vstack((y, y)).T
    sw_multi = np.vstack((sw, sw)).T
    score = _weighted_percentile(y_multi, sw_multi, 50)
    assert_array_equal(score, [1, 1])


def test_weighted_percentile_equal():
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert score == 0


def test_weighted_percentile_zero_weight():
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    score = _weighted_percentile(y, sw, 50)
    assert score == 1.0


def test_weighted_median_equal_weights():
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)

    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert_almost_equal(median, w_median)


def test_weighted_median_integer_weights():
    # Checks weighted percentile=0.5 is same as median when manually adjust
    # data for weights
    rng = np.random.RandomState(0)
    x = rng.randint(20, size=10)
    weights = rng.choice(5, size=10)
    x_manual = np.repeat(x, weights)

    median = np.median(x_manual)
    w_median = _weighted_percentile(x, weights)

    assert median == w_median
