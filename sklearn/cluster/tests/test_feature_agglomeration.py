"""
Tests for sklearn.cluster._feature_agglomeration
"""
# Authors: Sergul Aydore 2017
import numpy as np
from sklearn.cluster import FeatureAgglomeration
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_array_almost_equal


def test_feature_agglomeration():
    rng = np.random.RandomState(0)
    n_samples, n_features, n_clusters = 10000, 1000, 200
    X = rng.randn(n_samples, n_features)
    agglo_mean = FeatureAgglomeration(n_clusters=n_clusters,
                                      pooling_func=np.mean)
    agglo_median = FeatureAgglomeration(n_clusters=n_clusters,
                                        pooling_func=np.median)
    agglo_mean.fit(X)
    agglo_median.fit(X)
    assert_true(np.size(np.unique(agglo_mean.labels_)) == n_clusters)
    assert_true(np.size(np.unique(agglo_median.labels_)) == n_clusters)

    # Test transform
    X_red_mean = agglo_mean.transform(X)
    X_red_median = agglo_median.transform(X)
    assert_true(X_red_mean.shape[1] == n_clusters)
    assert_true(X_red_median.shape[1] == n_clusters)

    # Check that fitting with no samples raises a ValueError
    assert_raises(ValueError, agglo_mean.fit, X[:0])
    assert_raises(ValueError, agglo_median.fit, X[:0])

    # Test inverse transform
    X_full_mean = agglo_mean.inverse_transform(X_red_mean)
    X_full_median = agglo_mean.inverse_transform(X_red_median)
    assert_true(np.unique(X_full_mean[0]).size == n_clusters)
    assert_true(np.unique(X_full_median[0]).size == n_clusters)

    assert_array_almost_equal(agglo_mean.transform(X_full_mean),
                              X_red_mean)
    assert_array_almost_equal(agglo_mean.transform(X_full_median),
                              X_red_median)
