from warnings import simplefilter

import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_allclose

from sklearn.cluster import BisectingKMeans

import pytest


@pytest.mark.parametrize("bisecting_strategy", ["biggest_inertia", "largest_cluster"])
def test_three_clusters(bisecting_strategy):
    """Tries to perform bisect k-means for three clusters to check
    if splitting data is performed correctly
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
    """Test Bisecting K-Means with sparse data
    Checks if labels and centers are the same between dense and sparse
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
    assert_array_almost_equal(normal_centers, sparse_centers)


@pytest.mark.parametrize("n_clusters", [4, 5])
def test_n_clusters(n_clusters):
    """Test if resulting labels are in range [0, n_clusters - 1]"""

    rng = np.random.RandomState(0)
    X = rng.rand(10, 2)

    bisect_means = BisectingKMeans(n_clusters=n_clusters, random_state=0)
    bisect_means.fit(X)

    assert_array_equal(np.unique(bisect_means.labels_), np.arange(n_clusters))


def test_one_cluster():
    """Test warnings and performance for n_cluster = 1"""

    X = np.array([[1, 2], [10, 2], [10, 8]])

    with pytest.warns(None) as w:
        bisect_means = BisectingKMeans(n_clusters=1, random_state=0)
        bisect_means.fit(X)
        msg = (
            "BisectingKMeans might be inefficient for n_cluster smaller than 3 "
            + " - Use Normal KMeans from sklearn.cluster instead."
        )
        assert str(w[0].message) == msg
        print(bisect_means.labels_)

        # All labels from fit or predict should be equal 0
        assert all(bisect_means.predict(X) == 0)


@pytest.mark.parametrize(
    "param, match, single_value",
    [
        # Test bisecting_strategy param
        (
            {"bisecting_strategy": "None"},
            r"Bisect Strategy must be 'biggest_inertia', or 'largest_cluster'",
            False,
        ),
        # Test init array
        (
            {"init": np.ones((5, 2))},
            "BisectingKMeans does not support init as array.",
            False,
        ),
        # Test single X value
        (
            {"n_clusters": 1},
            # 'Found array with 0 sample(s) (shape=(0, 1))
            #  while a minimum of 1 is required by BisectingKMeans.'
            "a minimum of 1 is required by BisectingKMeans.",
            True,
        ),
    ],
)
def test_wrong_params(param, match, single_value):
    """Test Exceptions at check_params function"""

    simplefilter("ignore")

    if single_value:
        X = np.ones((0, 1))
    else:
        rng = np.random.RandomState(0)
        X = rng.rand(5, 2)

    with pytest.raises(ValueError, match=match):
        bisect_means = BisectingKMeans(n_clusters=3, n_init=1)
        bisect_means.set_params(**param)
        bisect_means.fit(X)


@pytest.mark.parametrize("is_sparse", [True, False])
def test_fit_predict(is_sparse):
    """Check if labels from fit(X) method are same as from fit(X).predict(X)"""
    rng = np.random.RandomState(0)

    X = rng.rand(10, 2)

    if is_sparse:
        X[X < 0.8] = 0
        X = sp.csr_matrix(X)

    bisect_means = BisectingKMeans(n_clusters=3, random_state=0)
    bisect_means.fit(X)

    assert_array_equal(bisect_means.labels_, bisect_means.predict(X))
