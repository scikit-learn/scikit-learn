import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import pytest

from sklearn.ensemble._hist_gradient_boosting.binning import (
    _BinMapper,
    _find_binning_thresholds as _find_binning_thresholds_orig,
    _map_to_bins
)
from sklearn.ensemble._hist_gradient_boosting.types import X_DTYPE
from sklearn.ensemble._hist_gradient_boosting.types import X_BINNED_DTYPE


DATA = np.random.RandomState(42).normal(
    loc=[0, 10], scale=[1, 0.01], size=(int(1e6), 2)
).astype(X_DTYPE)


def _find_binning_thresholds(data, max_bins=256, subsample=int(2e5),
                             random_state=None):
    # Just a redef to avoid having to pass arguments all the time (as the
    # function is private we don't use default values for parameters)
    return _find_binning_thresholds_orig(data, max_bins, subsample,
                                         random_state)


def test_find_binning_thresholds_regular_data():
    data = np.linspace(0, 10, 1001).reshape(-1, 1)
    bin_thresholds = _find_binning_thresholds(data, max_bins=10)
    assert_allclose(bin_thresholds[0], [1, 2, 3, 4, 5, 6, 7, 8, 9])
    assert len(bin_thresholds) == 1

    bin_thresholds = _find_binning_thresholds(data, max_bins=5)
    assert_allclose(bin_thresholds[0], [2, 4, 6, 8])
    assert len(bin_thresholds) == 1


def test_find_binning_thresholds_small_regular_data():
    data = np.linspace(0, 10, 11).reshape(-1, 1)

    bin_thresholds = _find_binning_thresholds(data, max_bins=5)
    assert_allclose(bin_thresholds[0], [2, 4, 6, 8])

    bin_thresholds = _find_binning_thresholds(data, max_bins=10)
    assert_allclose(bin_thresholds[0], [1, 2, 3, 4, 5, 6, 7, 8, 9])

    bin_thresholds = _find_binning_thresholds(data, max_bins=11)
    assert_allclose(bin_thresholds[0], np.arange(10) + .5)

    bin_thresholds = _find_binning_thresholds(data, max_bins=255)
    assert_allclose(bin_thresholds[0], np.arange(10) + .5)


def test_find_binning_thresholds_random_data():
    bin_thresholds = _find_binning_thresholds(DATA, random_state=0)
    assert len(bin_thresholds) == 2
    for i in range(len(bin_thresholds)):
        assert bin_thresholds[i].shape == (255,)  # 256 - 1
        assert bin_thresholds[i].dtype == DATA.dtype

    assert_allclose(bin_thresholds[0][[64, 128, 192]],
                    np.array([-0.7, 0.0, 0.7]), atol=1e-1)

    assert_allclose(bin_thresholds[1][[64, 128, 192]],
                    np.array([9.99, 10.00, 10.01]), atol=1e-2)


def test_find_binning_thresholds_low_n_bins():
    bin_thresholds = _find_binning_thresholds(DATA, max_bins=128,
                                              random_state=0)
    assert len(bin_thresholds) == 2
    for i in range(len(bin_thresholds)):
        assert bin_thresholds[i].shape == (127,)  # 128 - 1
        assert bin_thresholds[i].dtype == DATA.dtype


def test_find_binning_thresholds_invalid_n_bins():
    err_msg = 'no smaller than 2 and no larger than 256'
    with pytest.raises(ValueError, match=err_msg):
        _find_binning_thresholds(DATA, max_bins=1024)


def test_bin_mapper_n_features_transform():
    mapper = _BinMapper(max_bins=42, random_state=42).fit(DATA)
    err_msg = 'This estimator was fitted with 2 features but 4 got passed'
    with pytest.raises(ValueError, match=err_msg):
        mapper.transform(np.repeat(DATA, 2, axis=1))


@pytest.mark.parametrize('n_bins', [16, 128, 256])
def test_map_to_bins(n_bins):
    bin_thresholds = _find_binning_thresholds(DATA, max_bins=n_bins,
                                              random_state=0)
    binned = np.zeros_like(DATA, dtype=X_BINNED_DTYPE, order='F')
    _map_to_bins(DATA, bin_thresholds, binned)
    assert binned.shape == DATA.shape
    assert binned.dtype == np.uint8
    assert binned.flags.f_contiguous

    min_indices = DATA.argmin(axis=0)
    max_indices = DATA.argmax(axis=0)

    for feature_idx, min_idx in enumerate(min_indices):
        assert binned[min_idx, feature_idx] == 0
    for feature_idx, max_idx in enumerate(max_indices):
        assert binned[max_idx, feature_idx] == n_bins - 1


@pytest.mark.parametrize("n_bins", [5, 10, 42])
def test_bin_mapper_random_data(n_bins):
    n_samples, n_features = DATA.shape

    expected_count_per_bin = n_samples // n_bins
    tol = int(0.05 * expected_count_per_bin)

    mapper = _BinMapper(max_bins=n_bins, random_state=42).fit(DATA)
    binned = mapper.transform(DATA)

    assert binned.shape == (n_samples, n_features)
    assert binned.dtype == np.uint8
    assert_array_equal(binned.min(axis=0), np.array([0, 0]))
    assert_array_equal(binned.max(axis=0), np.array([n_bins - 1, n_bins - 1]))
    assert len(mapper.bin_thresholds_) == n_features
    for bin_thresholds_feature in mapper.bin_thresholds_:
        assert bin_thresholds_feature.shape == (n_bins - 1,)
        assert bin_thresholds_feature.dtype == DATA.dtype
    assert np.all(mapper.actual_n_bins_ == n_bins)

    # Check that the binned data is approximately balanced across bins.
    for feature_idx in range(n_features):
        for bin_idx in range(n_bins):
            count = (binned[:, feature_idx] == bin_idx).sum()
            assert abs(count - expected_count_per_bin) < tol


@pytest.mark.parametrize("n_samples, n_bins", [
    (5, 5),
    (5, 10),
    (5, 11),
    (42, 255)
])
def test_bin_mapper_small_random_data(n_samples, n_bins):
    data = np.random.RandomState(42).normal(size=n_samples).reshape(-1, 1)
    assert len(np.unique(data)) == n_samples

    mapper = _BinMapper(max_bins=n_bins, random_state=42)
    binned = mapper.fit_transform(data)

    assert binned.shape == data.shape
    assert binned.dtype == np.uint8
    assert_array_equal(binned.ravel()[np.argsort(data.ravel())],
                       np.arange(n_samples))


@pytest.mark.parametrize("n_bins, n_distinct, multiplier", [
    (5, 5, 1),
    (5, 5, 3),
    (255, 12, 42),
])
def test_bin_mapper_identity_repeated_values(n_bins, n_distinct, multiplier):
    data = np.array(list(range(n_distinct)) * multiplier).reshape(-1, 1)
    binned = _BinMapper(max_bins=n_bins).fit_transform(data)
    assert_array_equal(data, binned)


@pytest.mark.parametrize('n_distinct', [2, 7, 42])
def test_bin_mapper_repeated_values_invariance(n_distinct):
    rng = np.random.RandomState(42)
    distinct_values = rng.normal(size=n_distinct)
    assert len(np.unique(distinct_values)) == n_distinct

    repeated_indices = rng.randint(low=0, high=n_distinct, size=1000)
    data = distinct_values[repeated_indices]
    rng.shuffle(data)
    assert_array_equal(np.unique(data), np.sort(distinct_values))

    data = data.reshape(-1, 1)

    mapper_1 = _BinMapper(max_bins=n_distinct)
    binned_1 = mapper_1.fit_transform(data)
    assert_array_equal(np.unique(binned_1[:, 0]), np.arange(n_distinct))

    # Adding more bins to the mapper yields the same results (same thresholds)
    mapper_2 = _BinMapper(max_bins=min(256, n_distinct * 3))
    binned_2 = mapper_2.fit_transform(data)

    assert_allclose(mapper_1.bin_thresholds_[0], mapper_2.bin_thresholds_[0])
    assert_array_equal(binned_1, binned_2)


@pytest.mark.parametrize("n_bins, scale, offset", [
    (3, 2, -1),
    (42, 1, 0),
    (256, 0.3, 42),
])
def test_bin_mapper_identity_small(n_bins, scale, offset):
    data = np.arange(n_bins).reshape(-1, 1) * scale + offset
    binned = _BinMapper(max_bins=n_bins).fit_transform(data)
    assert_array_equal(binned, np.arange(n_bins).reshape(-1, 1))


@pytest.mark.parametrize('n_bins_small, n_bins_large', [
    (2, 2),
    (3, 3),
    (4, 4),
    (42, 42),
    (256, 256),
    (5, 17),
    (42, 256),
])
def test_bin_mapper_idempotence(n_bins_small, n_bins_large):
    assert n_bins_large >= n_bins_small
    data = np.random.RandomState(42).normal(size=30000).reshape(-1, 1)
    mapper_small = _BinMapper(max_bins=n_bins_small)
    mapper_large = _BinMapper(max_bins=n_bins_large)
    binned_small = mapper_small.fit_transform(data)
    binned_large = mapper_large.fit_transform(binned_small)
    assert_array_equal(binned_small, binned_large)


@pytest.mark.parametrize('max_bins', [10, 100, 256])
@pytest.mark.parametrize('diff', [-5, 0, 5])
def test_actual_n_bins(max_bins, diff):
    # Check that actual_n_bins is n_unique_values when
    # n_unique_values <= max_bins, else max_bins.

    n_unique_values = max_bins + diff
    X = list(range(n_unique_values)) * 2
    X = np.array(X).reshape(-1, 1)
    mapper = _BinMapper(max_bins=max_bins).fit(X)
    assert np.all(mapper.actual_n_bins_ == min(max_bins, n_unique_values))


def test_subsample():
    # Make sure bin thresholds are different when applying subsampling
    mapper_no_subsample = _BinMapper(subsample=None, random_state=0).fit(DATA)
    mapper_subsample = _BinMapper(subsample=256, random_state=0).fit(DATA)

    for feature in range(DATA.shape[1]):
        assert not np.allclose(mapper_no_subsample.bin_thresholds_[feature],
                               mapper_subsample.bin_thresholds_[feature],
                               rtol=1e-4)
