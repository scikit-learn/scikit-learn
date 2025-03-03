import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

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
    # Checks that `_weighted_percentile` and `np.median` (both at probability level=0.5
    # and with `sample_weights` being all 1s) return the same percentiles if the number
    # of the samples in the data is odd. In this special case, `_weighted_percentile`
    # always falls on a precise value (not on the next lower value) and is thus equal to
    # `np.median`.
    # As discussed in #17370, a similar check with an even number of samples does not
    # consistently hold, since then the lower of two percentiles might be selected,
    # while the median might lie in between.
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

    # Check when array and sample_weight boht 2D
    w2 = rng.choice(5, size=10)
    w_2d = np.vstack((w1, w2)).T

    w_median = _weighted_percentile(x_2d, w_2d)
    p_axis_0 = [
        _weighted_percentile(x_2d[:, i], w_2d[:, i]) for i in range(x_2d.shape[1])
    ]
    assert_allclose(w_median, p_axis_0)


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
