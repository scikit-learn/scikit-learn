"""Testing for K-means"""

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_raises

from ..k_means_ import KMeans, MiniBatchKMeans
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


def test_mini_batch_k_means_pp_init():
    np.random.seed(1)
    sample = X[0:X.shape[0] / 2]
    km = MiniBatchKMeans(init="random").partial_fit(sample)
    # Let's recalculate the inertia on the whole dataset
    km.partial_fit(X)
    inertia = km.inertia_
    km.partial_fit(X[X.shape[0] / 2:])
    # And again
    km.partial_fit(X)
    assert(km.inertia_ < inertia)


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

def test_mbk_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=init_array, n_init=1).fit(X)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]], k=n_clusters)
