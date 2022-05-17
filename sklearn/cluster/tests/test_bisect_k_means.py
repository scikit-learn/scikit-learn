import numpy as np
import pytest
import scipy.sparse as sp

from sklearn.utils._testing import assert_array_equal, assert_allclose
from sklearn.cluster import BisectingKMeans


@pytest.mark.parametrize("bisecting_strategy", ["biggest_inertia", "largest_cluster"])
def test_three_clusters(bisecting_strategy):
    """Tries to perform bisect k-means for three clusters to check
    if splitting data is performed correctly.
    """

    # X = np.array([[1, 2], [1, 4], [1, 0],
    #               [10, 2], [10, 4], [10, 0],
    #               [10, 6], [10, 8], [10, 10]])

    # X[0][1] swapped with X[1][1] intentionally for checking labeling
    X = np.array(
        [[1, 2], [10, 4], [1, 0], [10, 2], [1, 4], [10, 0], [10, 6], [10, 8], [10, 10]]
    )
    bisect_means = BisectingKMeans(
        n_clusters=3, random_state=0, bisecting_strategy=bisecting_strategy
    )
    bisect_means.fit(X)

    expected_centers = [[10, 2], [10, 8], [1, 2]]
    expected_predict = [2, 0]
    expected_labels = [2, 0, 2, 0, 2, 0, 1, 1, 1]

    assert_allclose(expected_centers, bisect_means.cluster_centers_)
    assert_array_equal(expected_predict, bisect_means.predict([[0, 0], [12, 3]]))
    assert_array_equal(expected_labels, bisect_means.labels_)


def test_sparse():
    """Test Bisecting K-Means with sparse data.

    Checks if labels and centers are the same between dense and sparse.
    """

    rng = np.random.RandomState(0)

    X = rng.rand(20, 2)
    X[X < 0.8] = 0
    X_csr = sp.csr_matrix(X)

    bisect_means = BisectingKMeans(n_clusters=3, random_state=0)

    bisect_means.fit(X_csr)
    sparse_centers = bisect_means.cluster_centers_

    bisect_means.fit(X)
    normal_centers = bisect_means.cluster_centers_

    # Check if results is the same for dense and sparse data
    assert_allclose(normal_centers, sparse_centers, atol=1e-8)


@pytest.mark.parametrize("n_clusters", [4, 5])
def test_n_clusters(n_clusters):
    """Test if resulting labels are in range [0, n_clusters - 1]."""

    rng = np.random.RandomState(0)
    X = rng.rand(10, 2)

    bisect_means = BisectingKMeans(n_clusters=n_clusters, random_state=0)
    bisect_means.fit(X)

    assert_array_equal(np.unique(bisect_means.labels_), np.arange(n_clusters))


def test_one_cluster():
    """Test single cluster."""

    X = np.array([[1, 2], [10, 2], [10, 8]])

    bisect_means = BisectingKMeans(n_clusters=1, random_state=0).fit(X)

    # All labels from fit or predict should be equal 0
    assert all(bisect_means.labels_ == 0)
    assert all(bisect_means.predict(X) == 0)

    assert_allclose(bisect_means.cluster_centers_, X.mean(axis=0).reshape(1, -1))


@pytest.mark.parametrize("is_sparse", [True, False])
def test_fit_predict(is_sparse):
    """Check if labels from fit(X) method are same as from fit(X).predict(X)."""
    rng = np.random.RandomState(0)

    X = rng.rand(10, 2)

    if is_sparse:
        X[X < 0.8] = 0
        X = sp.csr_matrix(X)

    bisect_means = BisectingKMeans(n_clusters=3, random_state=0)
    bisect_means.fit(X)

    assert_array_equal(bisect_means.labels_, bisect_means.predict(X))


@pytest.mark.parametrize("is_sparse", [True, False])
def test_dtype_preserved(is_sparse, global_dtype):
    """Check that centers dtype is the same as input data dtype."""
    rng = np.random.RandomState(0)
    X = rng.rand(10, 2).astype(global_dtype, copy=False)

    if is_sparse:
        X[X < 0.8] = 0
        X = sp.csr_matrix(X)

    km = BisectingKMeans(n_clusters=3, random_state=0)
    km.fit(X)

    assert km.cluster_centers_.dtype == global_dtype


@pytest.mark.parametrize("is_sparse", [True, False])
def test_float32_float64_equivalence(is_sparse):
    """Check that the results are the same between float32 and float64."""
    rng = np.random.RandomState(0)
    X = rng.rand(10, 2)

    if is_sparse:
        X[X < 0.8] = 0
        X = sp.csr_matrix(X)

    km64 = BisectingKMeans(n_clusters=3, random_state=0).fit(X)
    km32 = BisectingKMeans(n_clusters=3, random_state=0).fit(X.astype(np.float32))

    assert_allclose(km32.cluster_centers_, km64.cluster_centers_)
    assert_array_equal(km32.labels_, km64.labels_)
