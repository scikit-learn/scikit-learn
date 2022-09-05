"""
Tests for sklearn.cluster._feature_agglomeration
"""
# Authors: Sergul Aydore 2017
import numpy as np
import pytest
from numpy.testing import assert_array_equal
from sklearn.cluster import FeatureAgglomeration
from sklearn.utils._testing import assert_allclose, assert_array_almost_equal
from sklearn.datasets import make_blobs


def test_feature_agglomeration():
    n_clusters = 1
    X = np.array([0, 0, 1]).reshape(1, 3)  # (n_samples, n_features)

    agglo_mean = FeatureAgglomeration(n_clusters=n_clusters, pooling_func=np.mean)
    agglo_median = FeatureAgglomeration(n_clusters=n_clusters, pooling_func=np.median)
    agglo_mean.fit(X)
    agglo_median.fit(X)

    assert np.size(np.unique(agglo_mean.labels_)) == n_clusters
    assert np.size(np.unique(agglo_median.labels_)) == n_clusters
    assert np.size(agglo_mean.labels_) == X.shape[1]
    assert np.size(agglo_median.labels_) == X.shape[1]

    # Test transform
    Xt_mean = agglo_mean.transform(X)
    Xt_median = agglo_median.transform(X)
    assert Xt_mean.shape[1] == n_clusters
    assert Xt_median.shape[1] == n_clusters
    assert Xt_mean == np.array([1 / 3.0])
    assert Xt_median == np.array([0.0])

    # Test inverse transform
    X_full_mean = agglo_mean.inverse_transform(Xt_mean)
    X_full_median = agglo_median.inverse_transform(Xt_median)
    assert np.unique(X_full_mean[0]).size == n_clusters
    assert np.unique(X_full_median[0]).size == n_clusters

    assert_array_almost_equal(agglo_mean.transform(X_full_mean), Xt_mean)
    assert_array_almost_equal(agglo_median.transform(X_full_median), Xt_median)


def test_feature_agglomeration_feature_names_out():
    """Check `get_feature_names_out` for `FeatureAgglomeration`."""
    X, _ = make_blobs(n_features=6, random_state=0)
    agglo = FeatureAgglomeration(n_clusters=3)
    agglo.fit(X)
    n_clusters = agglo.n_clusters_

    names_out = agglo.get_feature_names_out()
    assert_array_equal(
        [f"featureagglomeration{i}" for i in range(n_clusters)], names_out
    )


def test_feature_agglomeration_numerical_consistency(global_random_seed):
    """Ensure numerical consistency among np.float32 and np.float64"""
    rng = np.random.RandomState(global_random_seed)
    X_64, _ = make_blobs(n_features=12, random_state=rng)
    X_32 = X_64.astype(np.float32)

    agglo_32 = FeatureAgglomeration(n_clusters=3)
    agglo_64 = FeatureAgglomeration(n_clusters=3)

    X_trans_64 = agglo_64.fit_transform(X_64)
    X_trans_32 = agglo_32.fit_transform(X_32)

    assert_allclose(X_trans_32, X_trans_64)
