import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from sklearn.utils._array_api import (
    yield_namespace_device_dtype_combinations,
    yield_namespaces,
)
from sklearn.utils.estimator_checks import _array_api_for_tests
from sklearn.utils.stats import _weighted_percentile


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)

    y = xp.empty(102, dtype=xp.float64, device=device)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = xp.ones(102, dtype=xp.float64, device=device)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 1


@pytest.mark.parametrize("array_namespace", yield_namespaces())
def test_weighted_percentile_equal(array_namespace):
    xp = _array_api_for_tests(array_namespace, device=None)
    y = xp.full(102, 0.0, dtype=xp.float64)
    sw = xp.ones(102, dtype=xp.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert score == 0


@pytest.mark.parametrize("array_namespace", yield_namespaces())
def test_weighted_percentile_zero_weight(array_namespace):
    xp = _array_api_for_tests(array_namespace, device=None)
    y = xp.full(102, 1.0, dtype=xp.float64)
    sw = xp.full(102, 0.0, dtype=xp.float64)
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 1.0


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_zero_weight_zero_percentile(
    array_namespace, device, dtype_name
):
    xp = _array_api_for_tests(array_namespace, device)
    y = xp.asarray([0, 1, 2, 3, 4, 5], device=device)
    sw = xp.asarray([0, 0, 1, 1, 1, 0], device=device)
    score = _weighted_percentile(y, sw, 0)
    assert approx(score) == 2

    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 3

    score = _weighted_percentile(y, sw, 100)
    assert approx(score) == 4


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_median_equal_weights(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    x = x.astype(dtype_name, copy=False)
    weights = xp.ones(x.shape)

    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert median == approx(w_median)


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_median_integer_weights(array_namespace, device, dtype_name):
    # Checks weighted percentile=0.5 is same as median when manually weight
    # data
    rng = np.random.RandomState(0)
    x = rng.randint(20, size=10)
    x = x.astype(dtype_name, copy=False)
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
