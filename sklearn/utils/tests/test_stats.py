import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.utils.stats import _weighted_percentile


def test_weighted_percentile(interpolation):
    y = np.empty(102, dtype=np.float64)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50, interpolation="nearest")
    assert score == pytest.approx(1)


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest"]
)
def test_weighted_percentile_equal(interpolation):
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50, interpolation=interpolation)
    assert score == 0


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest"]
)
def test_weighted_percentile_zero_weight(interpolation):
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    score = _weighted_percentile(y, sw, 50, interpolation=interpolation)
    assert pytest.approx(score) == 1.0


@pytest.mark.parametrize("interpolation", ["linear", "nearest"])
def test_weighted_median_equal_weights(interpolation):
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)

    median = np.median(x)
    w_median = _weighted_percentile(x, weights, interpolation=interpolation)
    assert median == pytest.approx(w_median)


@pytest.mark.parametrize("interpolation", ["linear", "nearest"])
def test_weighted_median_integer_weights(interpolation):
    # Checks weighted percentile=0.5 is same as median when manually weight
    # data
    rng = np.random.RandomState(0)
    x = rng.randint(20, size=10)
    weights = rng.choice(5, size=10)
    x_manual = np.repeat(x, weights)

    median = np.median(x_manual)
    w_median = _weighted_percentile(x, weights, interpolation=interpolation)

    assert median == pytest.approx(w_median)


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest"]
)
def test_weighted_percentile_2d(interpolation):
    # Check for when array 2D and sample_weight 1D
    rng = np.random.RandomState(0)
    x1 = rng.randint(10, size=10)
    w1 = rng.choice(5, size=10)

    x2 = rng.randint(20, size=10)
    x_2d = np.vstack((x1, x2)).T

    w_median = _weighted_percentile(x_2d, w1, interpolation=interpolation)
    p_axis_0 = [
        _weighted_percentile(x_2d[:, i], w1, interpolation=interpolation)
        for i in range(x_2d.shape[1])
    ]
    assert_allclose(w_median, p_axis_0)

    # Check when array and sample_weight both 2D
    w2 = rng.choice(5, size=10)
    w_2d = np.vstack((w1, w2)).T

    w_median = _weighted_percentile(x_2d, w_2d, interpolation=interpolation)
    p_axis_0 = [
        _weighted_percentile(x_2d[:, i], w_2d[:, i],
                             interpolation=interpolation)
        for i in range(x_2d.shape[1])
    ]
    assert_allclose(w_median, p_axis_0)


def test_weighted_percentile_np_median():
    # check that our weighted percentile lead to the same results than
    # unweighted NumPy implementation with unit weights for the median
    rng = np.random.RandomState(42)
    X = rng.randn(10)
    X.sort()
    sample_weight = np.ones(X.shape)

    np_median = np.median(X)
    sklearn_median = _weighted_percentile(
        X, sample_weight, percentile=50.0, interpolation="linear"
    )

    assert sklearn_median == pytest.approx(np_median)


@pytest.mark.parametrize("interpolation", ["linear", "lower", "higher"])
@pytest.mark.parametrize("percentile", np.arange(0, 101, 10))
def test_weighted_percentile_np_percentile(interpolation, percentile):
    rng = np.random.RandomState(0)
    X = rng.randn(10)
    X.sort()
    sample_weight = np.ones(X.shape) / X.shape[0]

    np_percentile = np.percentile(X, percentile, interpolation=interpolation)
    sklearn_percentile = _weighted_percentile(
        X, sample_weight, percentile=percentile, interpolation=interpolation,
    )

    assert sklearn_percentile == pytest.approx(np_percentile)


def test_weighted_percentile_wrong_interpolation():
    err_msg = "'interpolation' should be one of"
    with pytest.raises(ValueError, match=err_msg):
        X = np.random.randn(10)
        sample_weight = np.ones(X.shape)
        _weighted_percentile(X, sample_weight, 50, interpolation="xxxx")
