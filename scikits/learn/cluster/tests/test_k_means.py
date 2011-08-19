"""Testing for K-means"""

import numpy as np
from scipy import sparse as sp
from numpy.testing import assert_equal
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_raises
from nose.tools import assert_true

from ..k_means_ import KMeans, MiniBatchKMeans
from ...datasets.samples_generator import make_blobs

n_clusters = 3
centers = np.array([[1, 1], [-1, -1], [1, -1]]) + 10
X, _ = make_blobs(n_samples=60, n_features=2, centers=centers,
                  cluster_std=0.1, shuffle=False, random_state=0)

S = sp.csr_matrix(X)


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


def test_sparse_mini_batch_k_means_pp_init():
    np.random.seed(1)
    sample = S[0:S.shape[0] / 2]
    km = MiniBatchKMeans(init="random").partial_fit(sample)
    # Let's recalculate the inertia on the whole dataset
    km.partial_fit(S)
    inertia = km.inertia_
    km.partial_fit(S[S.shape[0] / 2:])
    # And again
    km.partial_fit(S)
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


def test_k_means_invalid_init():
    np.random.seed(1)
    k_means = KMeans(init="invalid", n_init=1)
    assert_raises(ValueError, k_means.fit, X, k=n_clusters)


def test_k_means_copyx():
    """Check if copy_x=False returns nearly equal X after de-centering."""
    np.random.seed(1)
    my_X = X.copy()
    k_means = KMeans(copy_x=False).fit(my_X, k=n_clusters)
    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)


def test_k_means_singleton():
    """Check k_means with bad initialization and singleton clustering."""
    np.random.seed(1)
    my_X = np.array([[1.1, 1.1], [0.9, 1.1], [1.1, 0.9], [0.9, 0.9]])
    array_init = np.array([[1.0, 1.0], [5.0, 5.0]])
    k_means = KMeans(init=array_init).fit(my_X, k=2)

    # must be singleton clustering
    assert_equal(np.unique(k_means.labels_).size, 1)


def test_mbk_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=init_array).fit(X)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]], k=n_clusters)


def test_sparse_mbk_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=init_array).fit(S)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]], k=n_clusters)


def test_sparse_mbk_means_pp_init():
    np.random.seed(1)
    mbk_means = MiniBatchKMeans(init="k-means++")
    assert_raises(ValueError, mbk_means.fit, S, k=n_clusters)


def test_sparse_mbk_means_callable_init():
    np.random.seed(1)

    def test_init(Xbar, k, random_state):
        return np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=test_init).fit(S)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]], k=n_clusters)


def test_k_means_fixed_array_init_fit():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = KMeans(init=init_array, n_init=1).fit(X, k=n_clusters)

    another_init_array = np.vstack([X[1], X[30], X[50]])
    other_k_means = KMeans(init=init_array, n_init=1)
    other_k_means.fit(X, init=another_init_array, k=n_clusters)
    assert_true(not np.allclose(k_means.init, other_k_means.init),
                "init attributes must be different")


def test_mbkm_fixed_array_init_fit():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = MiniBatchKMeans(init=init_array).fit(X, k=n_clusters)

    another_init_array = np.vstack([X[1], X[30], X[50]])
    other_k_means = MiniBatchKMeans(init=init_array)
    other_k_means.fit(X, init=another_init_array, k=n_clusters)
    assert_true(not np.allclose(k_means.init, other_k_means.init),
                "init attributes must be different")
