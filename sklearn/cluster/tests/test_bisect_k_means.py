import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal, assert_array_almost_equal

from sklearn.cluster import BisectKMeans

import pytest


@pytest.mark.parametrize("bisect_strategy", ["biggest_sse", "largest_cluster"])
def test_three_clusters(bisect_strategy):
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
    bisect_means = BisectKMeans(
        n_clusters=3, random_state=0, bisect_strategy=bisect_strategy
    )
    bisect_means.fit(X)

    expected_centers = [[1, 2], [10, 2], [10, 8]]
    expected_predict = [0, 1]
    expected_labels = [0, 1, 0, 1, 0, 1, 2, 2, 2]

    assert_array_equal(expected_centers, bisect_means.cluster_centers_)
    assert_array_equal(expected_predict, bisect_means.predict([[0, 0], [12, 3]]))
    assert_array_equal(expected_labels, bisect_means.labels_)


def test_sparse():
    """Test Bisecting K-Means with sparse data"""

    rng = np.random.RandomState(0)

    X = rng.rand(40, 2)
    X[X < 0.8] = 0
    X_csr = sp.csr_matrix(X)

    bisect_means = BisectKMeans(n_clusters=3, random_state=0)
    bisect_means.fit(X_csr)

    sparse_centers = bisect_means.cluster_centers_

    bisect_means.fit(X)
    normal_centers = bisect_means.cluster_centers_
    assert_array_almost_equal(normal_centers, sparse_centers)


# def test_perf():
#     """Test Bisecting K-Means with sparse data"""
#
#     rng = np.random.RandomState(0)
#     X = rng.rand(10, 2)
#     X = np.around(X, 2)
#
#     np.set_printoptions(precision=2)
#
#     bisect_means = BisectKMeans(n_clusters=3, random_state=0)
#     bisect_means.fit(X)
#     print(f"Labels:        {bisect_means.labels_}")
#     print(f"Predict:        {bisect_means.predict(X)}")
#     assert_array_almost_equal(bisect_means.labels_, bisect_means.predict(X))
