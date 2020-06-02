import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.utils.stats import _weighted_percentile


@pytest.mark.parametrize(
    "interpolation, expected_median",
    [("lower", 0), ("linear", 1), ("higher", 1)]
)
def test_weighted_percentile(interpolation, expected_median):
    y = np.empty(102, dtype=np.float64)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50, interpolation=interpolation)
    assert score == pytest.approx(expected_median)


@pytest.mark.parametrize("interpolation", ["linear", "lower", "higher"])
def test_weighted_percentile_equal(interpolation):
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50, interpolation=interpolation)
    assert score == 0


def test_weighted_median_equal_weights():
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)

    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert median == pytest.approx(w_median)


def test_weighted_median_integer_weights():
    # Checks weighted percentile=0.5 is same as median when manually weight
    # data
    rng = np.random.RandomState(0)
    x = rng.randint(20, size=10)
    weights = rng.choice(5, size=10)
    x_manual = np.repeat(x, weights)

    median = np.median(x_manual)
    w_median = _weighted_percentile(x, weights)

    assert median == pytest.approx(w_median)


@pytest.mark.parametrize("interpolation", ["linear", "lower", "higher"])
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
@pytest.mark.parametrize("percentile", np.arange(0, 101, 2.5))
def test_weighted_percentile_np_percentile(interpolation, percentile):
    rng = np.random.RandomState(0)
    X = rng.randn(10)
    X.sort()
    sample_weight = np.ones(X.shape)

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


@pytest.mark.parametrize("percentile", np.arange(2.5, 100, 2.5))
def test_weighted_percentile_non_unit_weight(percentile):
    # check the cumulative sum of the weight on the left and right side of the
    # percentile
    rng = np.random.RandomState(42)
    X = rng.randn(1000)
    X.sort()
    sample_weight = rng.random(X.shape)
    sample_weight = sample_weight / sample_weight.sum()
    sample_weight *= 100

    percentile_value = _weighted_percentile(X, sample_weight, percentile)
    X_percentile_idx = np.searchsorted(X, percentile_value)
    assert sample_weight[:X_percentile_idx - 1].sum() < percentile
    assert sample_weight[:X_percentile_idx + 1].sum() > percentile


@pytest.mark.parametrize("n_features", [None, 2])
@pytest.mark.parametrize("interpolation", ["linear", "higher", "lower"])
@pytest.mark.parametrize("percentile", np.arange(0, 101, 25))
def test_weighted_percentile_single_weight(n_features, interpolation,
                                           percentile):
    rng = np.random.RandomState(42)
    X = rng.randn(10) if n_features is None else rng.randn(10, n_features)
    X.sort(axis=0)
    sample_weight = np.zeros(X.shape)
    pos_weight_idx = 4
    sample_weight[pos_weight_idx] = 1

    percentile_value = _weighted_percentile(
        X, sample_weight, percentile=percentile, interpolation=interpolation
    )
    assert percentile_value == pytest.approx(X[pos_weight_idx])


@pytest.mark.parametrize("n_features", [None, 2])
def test_weighted_percentile_all_null_weight(n_features):
    rng = np.random.RandomState(42)
    X = rng.randn(10) if n_features is None else rng.randn(10, n_features)
    sample_weight = np.zeros(X.shape)

    err_msg = "All weights cannot be null when computing a weighted percentile"
    with pytest.raises(ValueError, match=err_msg):
        _weighted_percentile(X, sample_weight, 50)
