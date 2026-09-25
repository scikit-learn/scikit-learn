import math

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from sklearn.utils._mean_variance import _dense_mean_variance_axis0


def _reference(X):
    """Two-pass mean and variance with exactly rounded sums.

    `np.longdouble` is not used as it is just float64 on some platforms (e.g.
    Windows, macOS arm64).
    """
    X = X.astype(np.float64)
    n_samples = X.shape[0]
    mean = np.array([math.fsum(col) / n_samples for col in X.T])
    # For features with a large offset, `X - mean` is exact (Sterbenz lemma), and
    # the correction accounts for the rounding of `mean`.
    D = X - mean
    var = np.array(
        [math.fsum(d**2) / n_samples - (math.fsum(d) / n_samples) ** 2 for d in D.T]
    )
    return mean, var


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
# Shapes covering single/multiple tiles (4096 samples) and column strips (1024
# features), with remainders.
@pytest.mark.parametrize(
    "shape", [(1, 3), (7, 1), (100, 5), (10_000, 3), (50, 2500), (4096, 2)]
)
def test_dense_mean_variance_matches_numpy(dtype, order, shape, global_random_seed):
    rng = np.random.RandomState(global_random_seed)
    X = rng.normal(loc=3.0, scale=2.0, size=shape)
    X = np.asarray(X, dtype=dtype, order=order)

    mean, var = _dense_mean_variance_axis0(X)

    assert mean.dtype == var.dtype == np.float64
    assert mean.shape == var.shape == (shape[1],)
    X64 = X.astype(np.float64)
    assert_allclose(mean, X64.mean(axis=0), rtol=1e-12)
    assert_allclose(var, X64.var(axis=0), rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
def test_dense_mean_variance_constant_features(dtype, order):
    X = np.empty((10_000, 3), dtype=dtype, order=order)
    X[:, 0] = 0.0
    X[:, 1] = -1.5
    X[:, 2] = 1e6 + 0.1

    mean, var = _dense_mean_variance_axis0(X)

    assert_array_equal(mean, X[0].astype(np.float64))
    assert_array_equal(var, 0.0)


@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize(
    "dtype, offset",
    [(np.float32, 1e3), (np.float32, 1e5), (np.float64, 1e8), (np.float64, 1e10)],
)
def test_dense_mean_variance_large_offset(dtype, offset, order, global_random_seed):
    # Features with a large offset relative to their standard deviation: the
    # naive E[X²] - E[X]² formula loses all significant digits here.
    rng = np.random.RandomState(global_random_seed)
    X = offset + rng.normal(size=(20_000, 4))
    X = np.asarray(X, dtype=dtype, order=order)
    ref_mean, ref_var = _reference(X)

    mean, var = _dense_mean_variance_axis0(X)

    # The mean can't be more accurate than the spacing of floats around it.
    assert_allclose(mean, ref_mean, rtol=4 * np.finfo(np.float64).eps)
    assert_allclose(var, ref_var, rtol=1e-12)


@pytest.mark.parametrize("order", ["C", "F"])
def test_dense_mean_variance_adversarial_orderings(order, global_random_seed):
    # Sorted data (e.g. timestamps), and a first tile far away from the rest of
    # the data: the center used to accumulate a tile can be far from the tile
    # mean.
    rng = np.random.RandomState(global_random_seed)
    n_samples = 20_000
    sorted_feature = 1e9 + np.sort(rng.uniform(0, 1e3, size=n_samples))
    outlier_head = 1e8 + rng.normal(size=n_samples)
    outlier_head[:100] += 1e5
    X = np.asarray(np.stack([sorted_feature, outlier_head], axis=1), order=order)
    ref_mean, ref_var = _reference(X)

    mean, var = _dense_mean_variance_axis0(X)

    assert_allclose(mean, ref_mean, rtol=4 * np.finfo(np.float64).eps)
    assert_allclose(var, ref_var, rtol=1e-12)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_dense_mean_variance_memory_layouts(dtype, global_random_seed):
    rng = np.random.RandomState(global_random_seed)
    X = rng.normal(size=(300, 20)).astype(dtype)

    # Non-contiguous views are copied, C- and F-contiguous arrays are not.
    for X_view in [X, np.asfortranarray(X), X[::2], X[:, ::3], X[::2].T.copy().T]:
        mean, var = _dense_mean_variance_axis0(X_view)
        X64 = X_view.astype(np.float64)
        # Means are close to 0, and the standard deviation is 1.
        assert_allclose(mean, X64.mean(axis=0), atol=1e-14)
        assert_allclose(var, X64.var(axis=0), rtol=1e-12)


def test_dense_mean_variance_empty():
    mean, var = _dense_mean_variance_axis0(np.empty((0, 3)))
    assert_array_equal(mean, np.full(3, np.nan))
    assert_array_equal(var, np.full(3, np.nan))

    mean, var = _dense_mean_variance_axis0(np.empty((5, 0)))
    assert mean.shape == var.shape == (0,)


def test_dense_mean_variance_invalid_input():
    with pytest.raises(ValueError, match="Expected a 2D array"):
        _dense_mean_variance_axis0(np.ones(3))
    with pytest.raises(TypeError, match="Expected float32 or float64"):
        _dense_mean_variance_axis0(np.ones((3, 2), dtype=np.int64))
