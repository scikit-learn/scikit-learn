import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from sklearn._config import config_context
from sklearn.utils._array_api import (
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils.estimator_checks import _array_api_for_tests
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
    assert approx(score) == 1


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_array_api(array_namespace, device, dtype_name):
    with config_context(array_api_dispatch=True):
        xp = _array_api_for_tests(array_namespace, device)
        y_np = np.empty(102, dtype=dtype_name)
        y_np[:50] = 0
        y_np[-51:] = 2
        y_np[-1] = 100000
        y_np[50] = 1
        y = xp.asarray(y_np, device=device)
        sw_np = np.ones(102, dtype=dtype_name)
        sw = xp.asarray(sw_np, device=device)
        sw[-1] = 0.0
        score = _weighted_percentile(y, sw, 50)
        assert approx(score) == 1


def test_weighted_percentile_equal():
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 0


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_equal_array_api(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)
    y = xp.full(102, 0.0, dtype=xp.float64)
    sw = xp.ones(102, dtype=xp.float64)
    sw[-1] = 0.0
    with config_context(array_api_dispatch=True):
        score = _weighted_percentile(y, sw, 50)
        assert score == 0
"""


def test_weighted_percentile_zero_weight():
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 1.0


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_zero_weight_array_api(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)
    y = xp.full(102, 1.0, dtype=xp.float64)
    sw = xp.full(102, 0.0, dtype=xp.float64)
    with config_context(array_api_dispatch=True):
        score = _weighted_percentile(y, sw, 50)
        assert approx(score) == 1.0
"""


def test_weighted_percentile_zero_weight_zero_percentile():
    y = np.array([0, 1, 2, 3, 4, 5])
    sw = np.array([0, 0, 1, 1, 1, 0])
    score = _weighted_percentile(y, sw, 0)
    assert approx(score) == 2

    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 3

    score = _weighted_percentile(y, sw, 100)
    assert approx(score) == 4


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_zero_weight_zero_percentile_array_api(
    array_namespace, device, dtype_name
):
    xp = _array_api_for_tests(array_namespace, device)
    y = xp.asarray([0, 1, 2, 3, 4, 5], device=device)
    sw = xp.asarray([0, 0, 1, 1, 1, 0], device=device)

    with config_context(array_api_dispatch=True):
        score = _weighted_percentile(y, sw, 0)
        assert approx(score) == 2

        score = _weighted_percentile(y, sw, 50)
        assert approx(score) == 3

        score = _weighted_percentile(y, sw, 100)
        assert approx(score) == 4

"""


def test_weighted_median_equal_weights():
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)
    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert median == approx(w_median)


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_median_equal_weights_array_api(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)
    # Checks weighted percentile=0.5 is same as median when weights equal
    rng = np.random.RandomState(0)
    # Odd size as _weighted_percentile takes lower weighted percentile
    x = rng.randint(10, size=11)
    x = x.astype(dtype_name, copy=False)
    x = xp.asarray(x, device=device)
    weights = xp.ones(x.shape, device=device)
    median = np.median(x)

    with config_context(array_api_dispatch=True):
        w_median = _weighted_percentile(x, weights)
        assert median == approx(w_median)

"""


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


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_median_integer_weights_array_api(array_namespace, device, dtype_name):
    with config_context(array_api_dispatch=True):
        xp = _array_api_for_tests(array_namespace, device)
        # Checks weighted percentile=0.5 is same as median when manually weight
        # data
        rng = np.random.RandomState(0)
        x = rng.randint(20, size=10)
        x = x.astype(dtype_name, copy=False)
        x = xp.asarray(x, device=device)
        weights = rng.choice(5, size=10)
        weights = weights.astype(dtype_name, copy=False)
        weights = xp.asarray(weights, device=device)
        x_manual = xp.repeat(x, weights)

        median = np.median(x_manual)
        w_median = _weighted_percentile(x, weights)

        assert median == approx(w_median)

"""


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


"""
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_weighted_percentile_2d_array_api(array_namespace, device, dtype_name):
    with config_context(array_api_dispatch=True):
        xp = _array_api_for_tests(array_namespace, device)
        # Check for when array 2D and sample_weight 1D
        rng = np.random.RandomState(0)
        x1 = rng.randint(10, size=10)
        w1 = rng.choice(5, size=10)
        w1 = w1.astype(dtype_name, copy=False)
        w1 = xp.asarray(w1, device=device)

        x2 = rng.randint(20, size=10)
        x_2d = np.vstack((x1, x2)).T
        x_2d = x_2d.astype(dtype_name, copy=False)
        x_2d = xp.asarray(x_2d, device=device)

        w_median = _weighted_percentile(x_2d, w1)
        p_axis_0 = [_weighted_percentile(x_2d[:, i], w1) for i in range(x_2d.shape[1])]
        assert_allclose(w_median, p_axis_0)

        # Check when array and sample_weight boht 2D
        w2 = rng.choice(5, size=10)
        w_2d = np.vstack((w1, w2)).T
        w_2d = w_2d.astype(dtype_name, copy=False)
        w_2d = xp.asarray(w_2d, device=device)

        w_median = _weighted_percentile(x_2d, w_2d)
        p_axis_0 = [
            _weighted_percentile(x_2d[:, i], w_2d[:, i]) for i in range(x_2d.shape[1])
        ]
        assert_allclose(w_median, p_axis_0)
"""
