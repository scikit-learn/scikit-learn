"""
Tests for sklearn.cluster._feature_agglomeration
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.cluster import FeatureAgglomeration
from sklearn.datasets import make_blobs
from sklearn.utils._testing import assert_array_almost_equal


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


def test_feature_agglomeration_inverse_transform_array_like():
    """Check that inverse_transform works with lists, 1D arrays, and DataFrames.

    Non-regression test for issue #35052.
    """
    X = np.arange(20, dtype=float).reshape(10, 2)
    agglo = FeatureAgglomeration(n_clusters=1).fit(X)

    # 1. 2D list
    res_list_2d = agglo.inverse_transform([[0.0], [1.0]])
    assert isinstance(res_list_2d, np.ndarray)
    assert_array_equal(res_list_2d, [[0.0, 0.0], [1.0, 1.0]])

    # 2. 1D list
    res_list_1d = agglo.inverse_transform([0.5])
    assert isinstance(res_list_1d, np.ndarray)
    assert_array_equal(res_list_1d, [0.5, 0.5])

    # 3. 1D numpy array
    res_np_1d = agglo.inverse_transform(np.array([0.5]))
    assert isinstance(res_np_1d, np.ndarray)
    assert_array_equal(res_np_1d, [0.5, 0.5])

    # 4. pandas DataFrame and Series
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame([[0.0], [1.0]], columns=["c0"])
    res_df = agglo.inverse_transform(df)
    assert isinstance(res_df, np.ndarray)
    assert_array_equal(res_df, [[0.0, 0.0], [1.0, 1.0]])

    series = pd.Series([0.5])
    res_series = agglo.inverse_transform(series)
    assert isinstance(res_series, np.ndarray)
    assert_array_equal(res_series, [0.5, 0.5])

    # 5. set_output(transform="pandas") roundtrip
    agglo.set_output(transform="pandas")
    Xt_df = agglo.transform(X)
    assert isinstance(Xt_df, pd.DataFrame)
    res_set_output = agglo.inverse_transform(Xt_df)
    assert isinstance(res_set_output, np.ndarray)
    assert res_set_output.shape == (10, 2)
    assert_array_almost_equal(agglo.transform(res_set_output), Xt_df.to_numpy())


def test_feature_agglomeration_inverse_transform_shape_mismatch():
    """Check ValueError when feature count in inverse_transform does not match."""
    X = np.arange(20, dtype=float).reshape(10, 2)
    agglo = FeatureAgglomeration(n_clusters=1).fit(X)

    msg = "X has 2 features, but FeatureAgglomeration is expecting 1 features as input."
    with pytest.raises(ValueError, match=msg):
        agglo.inverse_transform([[1.0, 2.0]])

    with pytest.raises(ValueError, match=msg):
        agglo.inverse_transform([1.0, 2.0])


def test_feature_agglomeration_inverse_transform_dtype():
    """Check that inverse_transform preserves the input dtype."""
    X = np.arange(20, dtype=float).reshape(10, 2)
    agglo = FeatureAgglomeration(n_clusters=1).fit(X)

    X_int = np.array([[1], [2]], dtype=np.int32)
    res_int = agglo.inverse_transform(X_int)
    assert res_int.dtype == np.int32

    X_float32 = np.array([[1.0], [2.0]], dtype=np.float32)
    res_float32 = agglo.inverse_transform(X_float32)
    assert res_float32.dtype == np.float32
