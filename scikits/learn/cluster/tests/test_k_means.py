"""Testing for K-means"""

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_raises

from ..k_means_ import KMeans
from .common import generate_clustered_data

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters, std=.1)


def test_k_means_pp_init():
    np.random.seed(1)
    k_means = KMeans(init="k-means++").fit(X, k=n_clusters)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]], k=n_clusters)

def test_k_means_pp_random_init():
    np.random.seed(1)
    k_means = KMeans(init="random").fit(X, k=n_clusters)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]], k=n_clusters)


def test_k_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = KMeans(init=init_array, n_init=1).fit(X, k=n_clusters)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]], k=n_clusters)
