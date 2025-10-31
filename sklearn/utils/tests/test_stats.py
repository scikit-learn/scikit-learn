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
from sklearn.utils.stats import _weighted_percentile


@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("size", [10, 15])
def test_weighted_percentile_matches_median(size, average):
    """Ensure `_weighted_percentile` matches `median` when expected.

    With unit `sample_weight`, `_weighted_percentile` should match the median except
    when `average=False` and the number of samples is even.
    For an even array and `average=False`, `percentile_rank=50` gives the lower
    of the two 'middle' values, that are averaged when calculating the `median`.
    """
    y = np.arange(size)
    sample_weight = np.ones_like(y)

    score = _weighted_percentile(y, sample_weight, 50, average=average)

    # `_weighted_percentile(average=False)` does not match `median` when n is even
    if size % 2 == 0 and average is False:
        assert score != np.median(y)
    else:
        assert approx(score) == np.median(y)


@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("percentile_rank", [20, 35, 61, [5, 47]])
@pytest.mark.parametrize("size", [10, 15])
def test_weighted_percentile_matches_numpy(
    global_random_seed, size, percentile_rank, average
):
    """Check `_weighted_percentile` with unit weights is correct.

    `average=True` results should be the same as `np.percentile`'s
    'averaged_inverted_cdf'.
    `average=False` results should be the same as `np.percentile`'s
    'inverted_cdf'.
    Note `np.percentile` is the same as `np.quantile` except `q` is in range [0, 100].

    We parametrize through different `percentile_rank` and `size` to
    ensure we get cases where `g=0` and `g>0` (see Hyndman and Fan 1996 for details).
    """
    rng = np.random.RandomState(global_random_seed)
    y = rng.randint(20, size=size)
    sw = np.ones_like(y)

    score = _weighted_percentile(y, sw, percentile_rank, average=average)

    if average:
        method = "averaged_inverted_cdf"
    else:
        method = "inverted_cdf"

    assert approx(score) == np.percentile(y, percentile_rank, method=method)


@pytest.mark.parametrize("percentile_rank", [50, 100])
def test_weighted_percentile_plus_one_clip_max(percentile_rank):
    """Check `j+1` index is clipped to max, when `average=True`.

    `percentile_plus_one_indices` can exceed max index when `percentile_indices`
    is already at max index.
    Note that when `g` (Hyndman and Fan) / `fraction_above` is greater than 0,
    `j+1` (Hyndman and Fan) / `percentile_plus_one_indices` is calculated but
    never used, so it does not matter what this value is.
    When percentile of percentile rank 100 falls exactly on the last value in the
    `weighted_cdf`, `g=0` and `percentile_indices` is at max index. In this case
    we set `percentile_plus_one_indices` to be max index as well, so the result is
    the average of 2x the max index (i.e. last value of `weighted_cdf`).
    """
    # Note for both `percentile_rank`s 50 and 100,`percentile_indices` is already at
    # max index
    y = np.array([[0, 0], [1, 1]])
    sw = np.array([[0.1, 0.2], [2, 3]])
    score = _weighted_percentile(y, sw, percentile_rank, average=True)
    for idx in range(2):
        assert score[idx] == approx(1.0)


def test_weighted_percentile_equal():
    """Check `weighted_percentile` with unit weights and all 0 values in `array`."""
    y = np.zeros(102, dtype=np.float64)
    sw = np.ones(102, dtype=np.float64)
    score = _weighted_percentile(y, sw, 50)
    assert approx(score) == 0


# XXX: is this really what we want? Shouldn't we raise instead?
# https://github.com/scikit-learn/scikit-learn/issues/31032
def test_weighted_percentile_all_zero_weights():
    """Check `weighted_percentile` with all weights equal to 0 returns last index."""
    y = np.arange(10)
    sw = np.zeros(10)
    value = _weighted_percentile(y, sw, 50)
    assert approx(value) == 9.0


@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("percentile_rank, expected_value", [(0, 2), (50, 3), (100, 5)])
def test_weighted_percentile_ignores_zero_weight(
    average, percentile_rank, expected_value
):
    """Check leading, trailing and middle 0 weights behave correctly.

    Check that leading zero-weight observations are ignored when `percentile_rank=0`.
    See #20528 for details.
    Check that when `average=True` and the `j+1` ('plus one') index has sample weight
    of 0, it is ignored. Also check that trailing zero weight observations are ignored
    (e.g., when `percentile_rank=100`).
    """
    y = np.array([0, 1, 2, 3, 4, 5, 6])
    sw = np.array([0, 0, 1, 1, 0, 1, 0])

    value = _weighted_percentile(
        np.vstack((y, y)).T, np.vstack((sw, sw)).T, percentile_rank, average=average
    )
    for idx in range(2):
        assert approx(value[idx]) == expected_value


@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("percentile_rank", [20, 35, 50, 61])
def test_weighted_percentile_frequency_weight_semantics(
    global_random_seed, percentile_rank, average
):
    """Check integer weights give the same result as repeating values."""
    rng = np.random.RandomState(global_random_seed)
    x = rng.randint(20, size=10)
    weights = rng.choice(5, size=10)

    x_repeated = np.repeat(x, weights)
    percentile_weights = _weighted_percentile(
        x, weights, percentile_rank, average=average
    )
    percentile_repeated = _weighted_percentile(
        x_repeated, np.ones_like(x_repeated), percentile_rank, average=average
    )
    assert percentile_weights == approx(percentile_repeated)
    # Also check `percentile_rank=50` matches `median`
    if percentile_rank == 50 and average:
        assert percentile_weights == approx(np.median(x_repeated))


@pytest.mark.parametrize("constant", [5, 8])
@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("percentile_rank", [20, 35, 50, 61, [20, 35, 50, 61]])
def test_weighted_percentile_constant_multiplier(
    global_random_seed, percentile_rank, average, constant
):
    """Check multiplying weights by a constant does not change the result.

    Note scale invariance does not always hold when multiplying by a
    float due to cumulative sum numerical error (which grows proportional to n).
    """
    rng = np.random.RandomState(global_random_seed)
    x = rng.randint(20, size=20)
    weights = rng.choice(5, size=20)
    weights_multiplied = weights * constant

    percentile = _weighted_percentile(x, weights, percentile_rank, average=average)
    percentile_multiplier = _weighted_percentile(
        x, weights_multiplied, percentile_rank, average=average
    )
    assert percentile == approx(percentile_multiplier)


@pytest.mark.parametrize("percentile_rank", [50, [20, 35, 50]])
@pytest.mark.parametrize("average", [True, False])
def test_weighted_percentile_2d(global_random_seed, percentile_rank, average):
    """Check `_weighted_percentile` behaviour is correct when `array` is 2D."""
    # Check for when array 2D and sample_weight 1D
    rng = np.random.RandomState(global_random_seed)
    x1 = rng.randint(10, size=10)
    w1 = rng.choice(5, size=10)

    x2 = rng.randint(20, size=10)
    x_2d = np.vstack((x1, x2)).T

    wp = _weighted_percentile(
        x_2d, w1, percentile_rank=percentile_rank, average=average
    )

    if isinstance(percentile_rank, list):
        p_list = []
        for pr in percentile_rank:
            p_list.append(
                [
                    _weighted_percentile(
                        x_2d[:, i], w1, percentile_rank=pr, average=average
                    )
                    for i in range(x_2d.shape[1])
                ]
            )
        p_axis_0 = np.stack(p_list, axis=-1)
        assert wp.shape == (x_2d.shape[1], len(percentile_rank))
    else:
        # percentile_rank is scalar
        p_axis_0 = [
            _weighted_percentile(
                x_2d[:, i], w1, percentile_rank=percentile_rank, average=average
            )
            for i in range(x_2d.shape[1])
        ]
        assert wp.shape == (x_2d.shape[1],)

    assert_allclose(wp, p_axis_0)

    # Check when array and sample_weight both 2D
    w2 = rng.choice(5, size=10)
    w_2d = np.vstack((w1, w2)).T

    wp = _weighted_percentile(
        x_2d, w_2d, percentile_rank=percentile_rank, average=average
    )

    if isinstance(percentile_rank, list):
        p_list = []
        for pr in percentile_rank:
            p_list.append(
                [
                    _weighted_percentile(
                        x_2d[:, i], w_2d[:, i], percentile_rank=pr, average=average
                    )
                    for i in range(x_2d.shape[1])
                ]
            )
        p_axis_0 = np.stack(p_list, axis=-1)
        assert wp.shape == (x_2d.shape[1], len(percentile_rank))
    else:
        # percentile_rank is scalar
        p_axis_0 = [
            _weighted_percentile(
                x_2d[:, i], w_2d[:, i], percentile_rank=percentile_rank, average=average
            )
            for i in range(x_2d.shape[1])
        ]
        assert wp.shape == (x_2d.shape[1],)

    assert_allclose(wp, p_axis_0)


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize(
    "data, weights, percentile",
    [
        # NumPy scalars input (handled as 0D arrays on array API)
        (np.float32(42), np.int32(1), 50),
        # Random 1D array, constant weights
        (lambda rng: rng.rand(50), np.ones(50).astype(np.int32), 50),
        # Random 2D array and random 1D weights
        (lambda rng: rng.rand(50, 3), lambda rng: rng.rand(50).astype(np.float32), 75),
        # Random 2D array and random 2D weights
        (
            lambda rng: rng.rand(20, 3),
            lambda rng: rng.rand(20, 3).astype(np.float32),
            [25, 75],
        ),
        # zero-weights and `rank_percentile=0` (#20528) (`sample_weight` dtype: int64)
        (np.array([0, 1, 2, 3, 4, 5]), np.array([0, 0, 1, 1, 1, 0]), 0),
        # np.nan's in data and some zero-weights (`sample_weight` dtype: int64)
        (np.array([np.nan, np.nan, 0, 3, 4, 5]), np.array([0, 1, 1, 1, 1, 0]), 0),
        # `sample_weight` dtype: int32
        (
            np.array([0, 1, 2, 3, 4, 5]),
            np.array([0, 1, 1, 1, 1, 0], dtype=np.int32),
            [25, 75],
        ),
    ],
)
def test_weighted_percentile_array_api_consistency(
    global_random_seed, array_namespace, device, dtype_name, data, weights, percentile
):
    """Check `_weighted_percentile` gives consistent results with array API."""
    xp = _array_api_for_tests(array_namespace, device)

    # Skip test for percentile=0 edge case (#20528) on namespace/device where
    # xp.nextafter is broken. This is the case for torch with MPS device:
    # https://github.com/pytorch/pytorch/issues/150027
    zero = xp.zeros(1, device=device)
    one = xp.ones(1, device=device)
    if percentile == 0 and xp.all(xp.nextafter(zero, one) == zero):
        pytest.xfail(f"xp.nextafter is broken on {device}")

    rng = np.random.RandomState(global_random_seed)
    X_np = data(rng) if callable(data) else data
    weights_np = weights(rng) if callable(weights) else weights
    # Ensure `data` of correct dtype
    X_np = X_np.astype(dtype_name)

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

    # Check dtype correct (`sample_weight` should follow `array`)
    if dtype_name == "float32":
        assert result_xp_np.dtype == result_np.dtype == np.float32
    else:
        assert result_xp_np.dtype == np.float64


@pytest.mark.parametrize("average", [True, False])
@pytest.mark.parametrize("sample_weight_ndim", [1, 2])
def test_weighted_percentile_nan_filtered(
    global_random_seed, sample_weight_ndim, average
):
    """Test `_weighted_percentile` ignores NaNs.

    Calling `_weighted_percentile` on an array with nan values returns the same
    results as calling `_weighted_percentile` on a filtered version of the data.
    We test both with sample_weight of the same shape as the data and with
    one-dimensional sample_weight.
    """

    rng = np.random.RandomState(global_random_seed)
    array_with_nans = rng.rand(100, 10)
    array_with_nans[rng.rand(*array_with_nans.shape) < 0.5] = np.nan
    nan_mask = np.isnan(array_with_nans)

    if sample_weight_ndim == 2:
        sample_weight = rng.randint(1, 6, size=(100, 10))
    else:
        sample_weight = rng.randint(1, 6, size=(100,))

    # Find the weighted percentile on the array with nans:
    results = _weighted_percentile(array_with_nans, sample_weight, 30, average=average)

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
            _weighted_percentile(
                filtered_array[col], filtered_weights[col], 30, average=average
            )
            for col in range(array_with_nans.shape[1])
        ]
    )

    assert_array_equal(expected_results, results)


@pytest.mark.parametrize(
    "percentile_rank, expected",
    [
        (90, [np.nan, 5]),
        ([50, 90], [[np.nan, np.nan], [2.0, 5.0]]),
    ],
)
def test_weighted_percentile_all_nan_column(percentile_rank, expected):
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
    values = _weighted_percentile(array, weights, percentile_rank)

    # The percentile of the second column should be `5` even though there are many nan
    # values present; the percentile of the first column can only be nan, since there
    # are no other possible values:
    assert np.array_equal(values, expected, equal_nan=True)


@pytest.mark.skipif(
    np_version < parse_version("2.0"),
    reason="np.quantile only accepts weights since version 2.0",
)
@pytest.mark.parametrize("percentile", [66, 10, 50])
@pytest.mark.parametrize("average", [False, True])
@pytest.mark.parametrize("uniform_weight", [False, True])
def test_weighted_percentile_like_numpy_quantile(
    percentile, average, uniform_weight, global_random_seed
):
    """Check `_weighted_percentile` is equivalent to `np.quantile` with weights."""
    # TODO: remove the following skip once no longer applicable.
    if average and not uniform_weight:
        pytest.skip(
            "np.quantile does not support weights with method='averaged_inverted_cdf'"
        )

    rng = np.random.RandomState(global_random_seed)
    array = rng.rand(10, 100)
    if uniform_weight:
        sample_weight = np.ones_like(array) * rng.randint(1, 6, size=1)
    else:
        sample_weight = rng.randint(1, 6, size=(10, 100))

    percentile_weighted_percentile = _weighted_percentile(
        array, sample_weight, percentile, average=average
    )
    percentile_numpy_quantile = np.quantile(
        array,
        percentile / 100,
        weights=sample_weight if not uniform_weight else None,
        method="averaged_inverted_cdf" if average else "inverted_cdf",
        axis=0,
    )

    assert_array_equal(percentile_weighted_percentile, percentile_numpy_quantile)


@pytest.mark.skipif(
    np_version < parse_version("2.0"),
    reason="np.nanquantile only accepts weights since version 2.0",
)
@pytest.mark.parametrize("percentile", [66, 10, 50])
@pytest.mark.parametrize("average", [False, True])
@pytest.mark.parametrize("uniform_weight", [False, True])
def test_weighted_percentile_like_numpy_nanquantile(
    percentile, average, uniform_weight, global_random_seed
):
    """Check `_weighted_percentile` equivalent to `np.nanquantile` with weights."""
    # TODO: remove the following skip once no longer applicable.
    if average and not uniform_weight:
        pytest.skip(
            "np.nanquantile does not support weights with "
            "method='averaged_inverted_cdf'"
        )

    rng = np.random.RandomState(global_random_seed)
    array_with_nans = rng.rand(10, 100)
    array_with_nans[rng.rand(*array_with_nans.shape) < 0.5] = np.nan
    if uniform_weight:
        sample_weight = np.ones_like(array_with_nans) * rng.randint(
            1,
            6,
            size=1,
        )
    else:
        sample_weight = rng.randint(1, 6, size=(10, 100))

    percentile_weighted_percentile = _weighted_percentile(
        array_with_nans, sample_weight, percentile, average=average
    )
    percentile_numpy_nanquantile = np.nanquantile(
        array_with_nans,
        percentile / 100,
        weights=sample_weight if not uniform_weight else None,
        method="averaged_inverted_cdf" if average else "inverted_cdf",
        axis=0,
    )

    assert_array_equal(percentile_weighted_percentile, percentile_numpy_nanquantile)
