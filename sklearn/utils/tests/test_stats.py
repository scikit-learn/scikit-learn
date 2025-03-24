import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from sklearn._config import config_context
from sklearn.utils._array_api import (
    _convert_to_numpy,
    get_namespace,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._array_api import device as array_device
from sklearn.utils.estimator_checks import _array_api_for_tests
from sklearn.utils.fixes import np_version, parse_version
from sklearn.utils.stats import _averaged_weighted_percentile, _weighted_percentile


def test_averaged_weighted_median():
    y = np.array([0, 1, 2, 3, 4, 5])
    sw = np.array([1, 1, 1, 1, 1, 1])

    score = _averaged_weighted_percentile(y, sw, 50)

    assert score == np.median(y)


# TODO: remove @pytest.mark.skipif when numpy min version >= 1.22.
@pytest.mark.skipif(
    condition=np_version < parse_version("1.22"),
    reason="older numpy do not support the 'method' parameter",
)
def test_averaged_weighted_percentile():
    rng = np.random.RandomState(0)
    y = rng.randint(20, size=10)

    sw = np.ones(10)

    score = _averaged_weighted_percentile(y, sw, 20)

    assert score == np.percentile(y, 20, method="averaged_inverted_cdf")


def test_averaged_and_weighted_percentile():
    y = np.array([0, 1, 2])
    sw = np.array([5, 1, 5])
    q = 50

    score_averaged = _averaged_weighted_percentile(y, sw, q)
    score = _weighted_percentile(y, sw, q)

    assert score_averaged == score


def test_weighted_percentile():
    """Check `weighted_percentile` on artificial data with obvious median."""
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
    """Check `weighted_percentile` with all weights equal to 1."""
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 0


def test_weighted_percentile_zero_weight():
    """Check `weighted_percentile` with all weights equal to 0."""
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 1.0


def test_weighted_percentile_zero_weight_zero_percentile():
    """Check `weighted_percentile(percentile_rank=0)` behaves correctly.

    Ensures that (leading)zero-weight observations ignored when `percentile_rank=0`.
    See #20528 for details.
    """
    y = np.array([0, 1, 2, 3, 4, 5])
    sw = np.array([0, 0, 1, 1, 1, 0])
    value = _weighted_percentile(y, sw, 0)
    assert approx(value) == 2

    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 3

    value = _weighted_percentile(y, sw, 100)
    assert approx(value) == 4


def test_weighted_median_equal_weights():
    """Checks `_weighted_percentile(percentile_rank=50)` is the same as `np.median`.

    `sample_weights` are all 1s and the number of samples is odd.
    When number of samples is odd, `_weighted_percentile` always falls on a single
    observation (not between 2 values, in which case the lower value would be taken)
    and is thus equal to `np.median`.
    For an even number of samples, this check will not always hold as (note that
    for some other percentile methods it will always hold). See #17370 for details.
    """
    rng = np.random.RandomState(0)
    x = rng.randint(10, size=11)
    weights = np.ones(x.shape)
    median = np.median(x)
    w_median = _weighted_percentile(x, weights)
    assert median == approx(w_median)


def test_weighted_median_integer_weights():
    # Checks weighted percentile_rank=0.5 is same as median when manually weight
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
    # Check when array and sample_weight both 2D
    w2 = rng.choice(5, size=10)
    w_2d = np.vstack((w1, w2)).T

    w_median = _weighted_percentile(x_2d, w_2d)
    p_axis_0 = [
        _weighted_percentile(x_2d[:, i], w_2d[:, i]) for i in range(x_2d.shape[1])
    ]
    assert_allclose(w_median, p_axis_0)


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize(
    "data, weights, percentile",
    [
        # Random 1D array
        (lambda rng: rng.rand(50), np.ones(50), 50),
        # Random 1D array and random 1D weights
        (lambda rng: rng.rand(50), lambda rng: rng.rand(50), 75),
        # 2D array with 2D weights
        (lambda rng: rng.rand(20, 3), lambda rng: rng.rand(20, 3), 25),
        # zero-weights and `rank_percentile=0` (#20528)
        (np.array([0, 1, 2, 3, 4, 5]), np.array([0, 0, 1, 1, 1, 0]), 0),
    ],
)
def test_weighted_percentile_array_api_consistency(
    global_random_seed, array_namespace, device, dtype_name, data, weights, percentile
):
    """Check `_weighted_percentile` gives consistent results with array API."""
    if array_namespace == "array_api_strict":
        import array_api_strict

        if device == array_api_strict.Device("device1"):
            # See https://github.com/data-apis/array-api-strict/issues/134
            pytest.xfail(
                "array_api_strict bug with indexing with tuple of arrays "
                "on non-'CPU_DEVICE' devices."
            )

    xp = _array_api_for_tests(array_namespace, device)

    rng = np.random.RandomState(global_random_seed)
    X_np = data(rng) if callable(data) else data
    weights_np = weights(rng) if callable(weights) else weights
    # Ensure all inputs are the correct dtype
    X_np = X_np.astype(dtype_name)
    weights_np = weights_np.astype(dtype_name)

    result_np = _weighted_percentile(X_np, weights_np, percentile)
    # Convert to Array API arrays
    X_xp = xp.asarray(X_np, device=device)
    weights_xp = xp.asarray(weights_np, device=device)

    with config_context(array_api_dispatch=True):
        result_xp = _weighted_percentile(X_xp, weights_xp, percentile)
        assert array_device(result_xp) == array_device(X_xp)
        assert get_namespace(result_xp)[0] == get_namespace(X_xp)[0]
        result_xp_np = _convert_to_numpy(result_xp, xp=xp)

    assert result_xp_np.dtype == result_np.dtype
    assert result_xp_np.shape == result_np.shape
    assert_allclose(result_np, result_xp_np)


@pytest.mark.parametrize("sample_weight_ndim", [1, 2])
def test_weighted_percentile_nan_filtered(sample_weight_ndim):
    """Test that calling _weighted_percentile on an array with nan values returns
    the same results as calling _weighted_percentile on a filtered version of the data.
    We test both with sample_weight of the same shape as the data and with
    one-dimensional sample_weight."""

    rng = np.random.RandomState(42)
    array_with_nans = rng.rand(10, 100)
    array_with_nans[rng.rand(*array_with_nans.shape) < 0.5] = np.nan
    nan_mask = np.isnan(array_with_nans)

    if sample_weight_ndim == 2:
        sample_weight = rng.randint(1, 6, size=(10, 100))
    else:
        sample_weight = rng.randint(1, 6, size=(10,))

    # Find the weighted percentile on the array with nans:
    results = _weighted_percentile(array_with_nans, sample_weight, 30)

    # Find the weighted percentile on the filtered array:
    filtered_array = [
        array_with_nans[~nan_mask[:, col], col]
        for col in range(array_with_nans.shape[1])
    ]
    if sample_weight.ndim == 1:
        sample_weight = np.repeat(sample_weight, array_with_nans.shape[1]).reshape(
            array_with_nans.shape[0], array_with_nans.shape[1]
        )
    filtered_weights = [
        sample_weight[~nan_mask[:, col], col] for col in range(array_with_nans.shape[1])
    ]

    expected_results = np.array(
        [
            _weighted_percentile(filtered_array[col], filtered_weights[col], 30)
            for col in range(array_with_nans.shape[1])
        ]
    )

    assert_array_equal(expected_results, results)


def test_weighted_percentile_all_nan_column():
    """Check that nans are ignored in general, except for all NaN columns."""

    array = np.array(
        [
            [np.nan, 5],
            [np.nan, 1],
            [np.nan, np.nan],
            [np.nan, np.nan],
            [np.nan, 2],
            [np.nan, np.nan],
        ]
    )
    weights = np.ones_like(array)
    percentile_rank = 90

    values = _weighted_percentile(array, weights, percentile_rank)

    # The percentile of the second column should be `5` even though there are many nan
    # values present; the percentile of the first column can only be nan, since there
    # are no other possible values:
    assert np.array_equal(values, np.array([np.nan, 5]), equal_nan=True)


@pytest.mark.skipif(
    np_version < parse_version("2.0"),
    reason="np.quantile only accepts weights since version 2.0",
)
@pytest.mark.parametrize("percentile", [66, 10, 50])
def test_weighted_percentile_like_numpy_quantile(percentile):
    """Check that _weighted_percentile delivers equivalent results as np.quantile
    with weights."""

    rng = np.random.RandomState(42)
    array = rng.rand(10, 100)
    sample_weight = rng.randint(1, 6, size=(10, 100))

    percentile_weighted_percentile = _weighted_percentile(
        array, sample_weight, percentile
    )
    percentile_numpy_quantile = np.quantile(
        array, percentile / 100, weights=sample_weight, axis=0, method="inverted_cdf"
    )

    assert_array_equal(percentile_weighted_percentile, percentile_numpy_quantile)


@pytest.mark.skipif(
    np_version < parse_version("2.0"),
    reason="np.nanquantile only accepts weights since version 2.0",
)
@pytest.mark.parametrize("percentile", [66, 10, 50])
def test_weighted_percentile_like_numpy_nanquantile(percentile):
    """Check that _weighted_percentile delivers equivalent results as np.nanquantile
    with weights."""

    rng = np.random.RandomState(42)
    array_with_nans = rng.rand(10, 100)
    array_with_nans[rng.rand(*array_with_nans.shape) < 0.5] = np.nan
    sample_weight = rng.randint(1, 6, size=(10, 100))

    percentile_weighted_percentile = _weighted_percentile(
        array_with_nans, sample_weight, percentile
    )
    percentile_numpy_nanquantile = np.nanquantile(
        array_with_nans,
        percentile / 100,
        weights=sample_weight,
        axis=0,
        method="inverted_cdf",
    )

    assert_array_equal(percentile_weighted_percentile, percentile_numpy_nanquantile)
