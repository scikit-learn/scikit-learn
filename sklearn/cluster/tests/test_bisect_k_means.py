import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal, assert_array_almost_equal

from sklearn.cluster import BisectKMeans

import pytest


@pytest.mark.parametrize("bisect_strategy", ["biggest_sse",
                                             "child_biggest_sse",
                                             "largest_cluster"])
def test_three_clusters(bisect_strategy):
    """ Tries to perform bisect k-means for three clusters to check
        if splitting data is performed correctly
    """

    # X = np.array([[1, 2], [1, 4], [1, 0],
    #               [10, 2], [10, 4], [10, 0],
    #               [10, 6], [10, 8], [10, 10]])

    # X[0][1] swapped with X[1][1] intentionally for checking labeling

    X = np.array([[1, 2], [10, 4], [1, 0],
                  [10, 2], [1, 4], [10, 0],
                  [10, 6], [10, 8], [10, 10]])
    bisect_means = BisectKMeans(n_clusters=3, random_state=0,
                                bisect_strategy=bisect_strategy)
    bisect_means.fit(X)

    # "largest_cluster" should produce same results
    # but with different ordering
    if bisect_strategy == "largest_cluster":
        expected_centers = [[1, 2], [10, 2], [10, 8]]
        expected_predict = [0, 1]
        expected_labels = [0, 1, 0, 1, 0, 1, 2, 2, 2]
    else:
        expected_centers = [[1, 2], [10, 8], [10, 2]]
        expected_predict = [0, 2]
        expected_labels = [0, 2, 0, 2, 0, 2, 1, 1, 1]

    assert_array_equal(expected_centers, bisect_means.cluster_centers_)
    assert_array_equal(bisect_means.predict([[0, 0], [12, 3]]),
                       expected_predict)
    assert_array_equal(bisect_means.labels_, expected_labels)


def test_sparse():
    """ Test Bisecting K-Means with sparse data """

    rng = np.random.RandomState(0)

    X = rng.rand(40, 2)
    X[X < .8] = 0
    X_csr = sp.csr_matrix(X)

    bisect_means = BisectKMeans(n_clusters=3, random_state=0)
    bisect_means.fit(X_csr)

    sparse_centers = bisect_means.cluster_centers_

    bisect_means.fit(X)
    normal_centers = bisect_means.cluster_centers_
    assert_array_almost_equal(normal_centers, sparse_centers)


@pytest.mark.filterwarnings("ignore:Explicit initial center*:RuntimeWarning")
def test_init_array():
    """ Test Bisecting K-Means with init array
    Note that it would work only for bisect_strategy= 'child_biggest_sse'
    Other strategies doesn't support init as array
    """

    X = np.array([[1, 2], [1, 4], [1, 0],
                  [10, 2], [10, 4], [10, 0],
                  [10, 6], [10, 8], [10, 10]])
    init = np.array([[1, 1], [5, 6], [10, 2]])
    bisect_means = BisectKMeans(n_clusters=3, random_state=0, init=init,
                                bisect_strategy="child_biggest_sse")
    bisect_means.fit(X)