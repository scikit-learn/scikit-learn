import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from sklearn.utils.stats import _weighted_percentile


def test_weighted_percentile():
    y = np.empty(102, dtype=np.float64)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 1


def test_weighted_percentile_equal():
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    value = _weighted_percentile(y, sw, 50)
    assert value == 0


def test_weighted_percentile_zero_weight():
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 1.0


def test_weighted_percentile_zero_weight_zero_percentile():
    y = np.array([0, 1, 2, 3, 4, 5])
    sw = np.array([0, 0, 1, 1, 1, 0])
    value = _weighted_percentile(y, sw, 0)
    assert approx(value) == 2

    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 3

    value = _weighted_percentile(y, sw, 100)
    assert approx(value) == 4


def test_weighted_median_equal_weights():
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)

    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert median == approx(w_median)


def test_weighted_median_integer_weights():
    # Checks weighted percentile=0.5 is same as median when manually weight
    # data
    rng = np.random.RandomState(0)
    x = rng.randint(20, size=10)
    weights = rng.choice(5, size=10)
    x_manual = np.repeat(x, weights)

    median = np.median(x_manual)
    w_median = _weighted_percentile(x, weights)

    assert median == approx(w_median)


def test_weighted_percentile_2d():
    # Check for when array 2D and sample_weight 1D
    rng = np.random.RandomState(0)
    x1 = rng.randint(10, size=10)
    w1 = rng.choice(5, size=10)

    x2 = rng.randint(20, size=10)
    x_2d = np.vstack((x1, x2)).T

    w_median = _weighted_percentile(x_2d, w1)
    p_axis_0 = [_weighted_percentile(x_2d[:, i], w1) for i in range(x_2d.shape[1])]
    assert_allclose(w_median, p_axis_0)

    # Check when array and sample_weight boht 2D
    w2 = rng.choice(5, size=10)
    w_2d = np.vstack((w1, w2)).T

    w_median = _weighted_percentile(x_2d, w_2d)
    p_axis_0 = [
        _weighted_percentile(x_2d[:, i], w_2d[:, i]) for i in range(x_2d.shape[1])
    ]
    assert_allclose(w_median, p_axis_0)


def test_weighted_percentile_nan():
    """Test that calling _weighted_percentile on an array with nan values returns
    the same results as calling _weighted_percentile on a filtered version of the data.
    We test both with sample_weight of the same shape as the data and for
    one-dimensional sample_weight."""

    rng = np.random.RandomState(42)
    array = rng.rand(10, 100)

    nan_array = array.copy()
    nan_array[rng.rand(*nan_array.shape) < 0.5] = np.nan
    nan_mask = np.isnan(nan_array)

    sample_weights = [rng.randint(1, 6, size=(10, 100)), rng.randint(1, 6, size=(10,))]

    for sample_weight in sample_weights:
        # Calculate the weighted percentile on the array with nans:
        values_nan = _weighted_percentile(nan_array, sample_weight, 30)

        # Calculate the weighted percentile on the filtered array:
        filtered_values = [
            array[~nan_mask[:, col], col] for col in range(array.shape[1])
        ]
        if sample_weight.ndim == 1:
            sample_weight = np.repeat(sample_weight, array.shape[1]).reshape(
                array.shape[0], array.shape[1]
            )
        filtered_weights = [
            sample_weight[~nan_mask[:, col], col] for col in range(array.shape[1])
        ]

        values_filtered = np.array(
            [
                _weighted_percentile(filtered_values[col], filtered_weights[col], 30)
                for col in range(array.shape[1])
            ]
        )

        assert_array_equal(values_filtered, values_nan)
