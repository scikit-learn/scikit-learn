"""
Tests for sklearn.cluster._feature_agglomeration
"""
# Authors: Sergul Aydore 2017
import numpy as np
from sklearn.cluster import FeatureAgglomeration
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_almost_equal


def test_feature_agglomeration():
    n_clusters = 1
    X = np.array([0, 0, 1], ndmin=2)  # (n_samples, n_features)

    agglo_mean = FeatureAgglomeration(n_clusters=n_clusters,
                                      pooling_func=np.mean)
    agglo_median = FeatureAgglomeration(n_clusters=n_clusters,
                                        pooling_func=np.median)
    agglo_mean.fit(X)
    agglo_median.fit(X)
    assert_true(np.size(np.unique(agglo_mean.labels_)) == n_clusters)
    assert_true(np.size(np.unique(agglo_median.labels_)) == n_clusters)

    # Test transform
    Xt_mean = agglo_mean.transform(X)
    Xt_median = agglo_median.transform(X)
    assert_true(Xt_mean.shape[1] == n_clusters)
    assert_true(Xt_median.shape[1] == n_clusters)

    # Test inverse transform
    X_full_mean = agglo_mean.inverse_transform(Xt_mean)
    X_full_median = agglo_median.inverse_transform(Xt_median)
    assert_true(np.unique(X_full_mean[0]).size == n_clusters)
    assert_true(np.unique(X_full_median[0]).size == n_clusters)

    assert_array_almost_equal(agglo_mean.transform(X_full_mean),
                              Xt_mean)
    assert_array_almost_equal(agglo_median.transform(X_full_median),
                              Xt_median)
