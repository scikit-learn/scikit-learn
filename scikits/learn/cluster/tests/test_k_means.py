"""Testing for K-means"""

import numpy as np
from scipy import sparse as sp
from numpy.testing import assert_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_raises
from nose.tools import assert_true

from ..k_means_ import KMeans, MiniBatchKMeans
from ...datasets.samples_generator import make_blobs
from .common import generate_clustered_data
from ...utils import shuffle

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters, std=.1)
X_csr = sp.csr_matrix(X)


def test_k_means_pp_init():
    np.random.seed(1)
    k_means = KMeans(init="k-means++", k=n_clusters).fit(X)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]])


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
    sample = X_csr[0:X_csr.shape[0] / 2]
    km = MiniBatchKMeans(init="random").partial_fit(sample)
    # Let's recalculate the inertia on the whole dataset
    km.partial_fit(X_csr)
    inertia = km.inertia_
    km.partial_fit(X_csr[X_csr.shape[0] / 2:])
    # And again
    km.partial_fit(X_csr)
    assert(km.inertia_ < inertia)


def test_k_means_pp_random_init():
    np.random.seed(1)
    k_means = KMeans(init="random", k=n_clusters).fit(X)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]])


def test_k_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = KMeans(init=init_array, n_init=1, k=n_clusters).fit(X)

    centers = k_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = k_means.labels_
    assert_equal(np.unique(labels).size, 3)
    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)

    # check error on dataset being too small
    assert_raises(ValueError, k_means.fit, [[0., 1.]])


def test_k_means_invalid_init():
    np.random.seed(1)
    k_means = KMeans(init="invalid", n_init=1, k=n_clusters)
    assert_raises(ValueError, k_means.fit, X)


def test_k_means_copyx():
    """Check if copy_x=False returns nearly equal X after de-centering."""
    np.random.seed(1)
    my_X = X.copy()
    k_means = KMeans(copy_x=False, k=n_clusters).fit(my_X)
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
    k_means = KMeans(init=array_init, k=2).fit(my_X)

    # must be singleton clustering
    assert_equal(np.unique(k_means.labels_).size, 1)


def test_mbk_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=init_array, k=n_clusters).fit(X)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]])


def test_sparse_mbk_means_fixed_array_init():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=init_array).fit(X_csr)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]])


def test_sparse_mbk_means_pp_init():
    np.random.seed(1)
    mbk_means = MiniBatchKMeans(init="k-means++", k=n_clusters)
    assert_raises(ValueError, mbk_means.fit, X_csr)


def test_sparse_mbk_means_callable_init():
    np.random.seed(1)

    def test_init(Xbar, k, random_state):
        return np.vstack([X[5], X[25], X[45]])
    mbk_means = MiniBatchKMeans(init=test_init).fit(X_csr)

    centers = mbk_means.cluster_centers_
    assert_equal(centers.shape, (n_clusters, 2))

    labels = mbk_means.labels_
    assert_equal(np.unique(labels).size, 3)

    assert_raises(ValueError, mbk_means.fit, [[0., 1.]])


def test_k_means_fixed_array_init_fit():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = KMeans(init=init_array, n_init=1, k=n_clusters).fit(X)

    another_init_array = np.vstack([X[1], X[30], X[50]])
    other_k_means = KMeans(init=init_array, n_init=1, k=n_clusters)
    other_k_means.set_params(init=another_init_array)
    other_k_means.fit(X)
    assert_true(not np.allclose(k_means.init, other_k_means.init),
                "init attributes must be different")


def test_mbkm_fixed_array_init_fit():
    np.random.seed(1)
    init_array = np.vstack([X[5], X[25], X[45]])
    k_means = MiniBatchKMeans(init=init_array, k=n_clusters).fit(X)

    another_init_array = np.vstack([X[1], X[30], X[50]])
    other_k_means = MiniBatchKMeans(init=another_init_array, k=n_clusters)
    other_k_means.fit(X)
    assert_true(not np.allclose(k_means.init, other_k_means.init),
                "init attributes must be different")


def test_mbk_means():
    n_samples = X.shape[0]
    true_labels = np.zeros((n_samples,), dtype=np.int32)
    true_labels[20:40] = 1
    true_labels[40:] = 2
    chunk_size = n_samples / 10
    # make sure init clusters are in different clusters
    init_array = np.vstack([X[5], X[25], X[45]])

    # shuffle original data
    X_shuffled, X_csr_shuffled, true_labels = shuffle(X, X_csr, true_labels,
                                                      random_state=1)

    mbk_means = MiniBatchKMeans(init=init_array, chunk_size=chunk_size,
                                k=n_clusters, random_state=1)
    mbk_means.fit(X_shuffled)
    assert_equal(true_labels, mbk_means.labels_)

    mbk_means = MiniBatchKMeans(init=init_array, chunk_size=chunk_size,
                                k=n_clusters, random_state=1)
    mbk_means.fit(X_csr_shuffled)
    assert_equal(true_labels, mbk_means.labels_)


def test_predict():
    k_means = KMeans(k=n_clusters).fit(X)

    # sanity check: predict centroid labels
    pred = k_means.predict(k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = k_means.predict(X)
    assert_array_equal(k_means.predict(X), k_means.labels_)


def test_predict_minibatch():
    mbk_means = MiniBatchKMeans(k=n_clusters).fit(X)

    # sanity check: predict centroid labels
    pred = mbk_means.predict(mbk_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = mbk_means.predict(X)
    assert_array_equal(mbk_means.predict(X), mbk_means.labels_)


def test_predict_minibatch_sparse_input():
    mbk_means = MiniBatchKMeans(k=n_clusters).fit(X_csr)

    # sanity check: predict centroid labels
    pred = mbk_means.predict(mbk_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = mbk_means.predict(X_csr)
    assert_array_equal(mbk_means.predict(X), mbk_means.labels_)

    # check that models trained on sparse input also works for dense input at
    # predict time
    pred = mbk_means.predict(X)
    assert_array_equal(mbk_means.predict(X), mbk_means.labels_)


def test_transform():
    k_means = KMeans(k=n_clusters)
    k_means.fit(X)
    X_new = k_means.transform(k_means.cluster_centers_)

    for c in range(n_clusters):
        assert_equal(X_new[c, c], 0)
        for c2 in range(n_clusters):
            if c != c2:
                assert_true(X_new[c, c2] > 0)
