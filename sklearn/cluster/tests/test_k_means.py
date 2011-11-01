"""Testing for K-means"""

import numpy as np
from scipy import sparse as sp
from numpy.testing import assert_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import assert_almost_equal

from sklearn.metrics.cluster import v_measure_score
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster.k_means_ import _labels_inertia
from sklearn.cluster._k_means import _mini_batch_update_csr
from sklearn.cluster._k_means import csr_row_norm_l2
from sklearn.datasets.samples_generator import make_blobs


# non centered, sparse centers to check the
centers = np.array([
    [0.0, 5.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 4.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 5.0, 1.0],
])
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)

X_csr = sp.csr_matrix(X)


def test_square_norms():
    x_squared_norms = (X ** 2).sum(axis=1)
    x_squared_norms_from_csr = csr_row_norm_l2(X_csr)
    assert_array_almost_equal(x_squared_norms,
                              x_squared_norms_from_csr, 5)


def test_labels_assignement_and_inertia():
    # pure numpy implementation as easily auditable reference gold
    # implementation
    rng = np.random.RandomState(42)
    noisy_centers = centers + rng.normal(size=centers.shape)
    labels_gold = - np.ones(n_samples, dtype=np.int)
    mindist = np.empty(n_samples)
    mindist.fill(np.infty)
    for center_id in range(n_clusters):
        dist = np.sum((X - centers[center_id]) ** 2, axis=1)
        labels_gold[dist < mindist] = center_id
        mindist = np.minimum(dist, mindist)
    inertia_gold = mindist.sum()
    assert_true((labels_gold != -1).all())

    # perform label assignement using the dense array input
    x_squared_norms = (X ** 2).sum(axis=1)
    labels_array, inertia_array = _labels_inertia(
        X, x_squared_norms, noisy_centers)
    assert_array_almost_equal(inertia_array, inertia_gold, 7)
    assert_array_equal(labels_array, labels_gold)

    # perform label assignement using the sparse CSR input
    x_squared_norms_from_csr = csr_row_norm_l2(X_csr)
    labels_csr, inertia_csr = _labels_inertia(
        X_csr, x_squared_norms_from_csr, noisy_centers)
    assert_array_almost_equal(inertia_csr, inertia_gold, 7)
    assert_array_equal(labels_csr, labels_gold)


def test_minibatch_update_sparse_csr():
    rng = np.random.RandomState(42)
    old_centers = centers + rng.normal(size=centers.shape)
    new_centers = old_centers.copy()
    counts = np.ones(new_centers.shape[0], dtype=np.int32)
    x_squared_norms = csr_row_norm_l2(X_csr, squared=True)
    buffer = np.zeros(centers.shape[1], dtype=np.double)

    slice_ = slice(0, 10)
    old_inertia, incremental_diff = _mini_batch_update_csr(
        X_csr, x_squared_norms, slice_, new_centers, counts, buffer, 1)
    assert_true(old_inertia > 0.0)

    # compute the new inertia on the same batch to check that is decreased
    _, new_inertia = _labels_inertia(X_csr[slice_], x_squared_norms,
                                     new_centers)
    assert_true(new_inertia > 0.0)
    assert_true(new_inertia < old_inertia)

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers - old_centers) ** 2)
    assert_almost_equal(incremental_diff, effective_diff, 5)


def _check_fitted_model(km):
    centers = km.cluster_centers_
    assert_equal(centers.shape, (n_clusters, n_features))

    labels = km.labels_
    assert_equal(np.unique(labels).shape[0], n_clusters)

    # check that the labels assignements are perfect (up to a permutation)
    assert_equal(v_measure_score(true_labels, labels), 1.0)
    assert_true(km.inertia_ > 0.0)

    # check error on dataset being too small
    assert_raises(ValueError, km.fit, [[0., 1.]])


def test_k_means_plus_plus_init():
    k_means = KMeans(init="k-means++", k=n_clusters, random_state=42).fit(X)
    _check_fitted_model(k_means)


def test_k_means_random_init():
    k_means = KMeans(init="random", k=n_clusters, random_state=42).fit(X)
    _check_fitted_model(k_means)


def test_k_means_perfect_init():
    k_means = KMeans(init=centers.copy(), k=n_clusters, random_state=42,
                     n_init=1)
    k_means.fit(X)
    _check_fitted_model(k_means)


def test_mb_k_means_plus_plus_init_dense_array():
    mb_k_means = MiniBatchKMeans(init="k-means++", k=n_clusters,
                                random_state=42)
    mb_k_means.fit(X)
    _check_fitted_model(mb_k_means)


def test_mb_k_means_plus_plus_init_sparse_csr_invalid():
    mb_k_means = MiniBatchKMeans(init="k-means++", k=n_clusters)
    assert_raises(ValueError, mb_k_means.fit, X_csr)


def test_minibatch_k_means_random_init_dense_array():
    mb_k_means = MiniBatchKMeans(init="random", k=n_clusters,
                                 random_state=42).fit(X)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_random_init_sparse_csr():
    mb_k_means = MiniBatchKMeans(init="random", k=n_clusters,
                                 random_state=42).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_perfect_init_dense_array():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), k=n_clusters,
                                 random_state=42).fit(X)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_perfect_init_sparse_csr():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), k=n_clusters,
                                 random_state=42).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_sparse_mb_k_means_callable_init():

    def test_init(X, k, random_state):
        return centers

    mb_k_means = MiniBatchKMeans(init=test_init, random_state=42).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_mini_batch_k_means_random_init_partial_fit():
    km = MiniBatchKMeans(k=n_clusters, init="random", random_state=42)

    # use the partial_fit API for online learning
    for X_minibatch in np.array_split(X, 10):
        km.partial_fit(X_minibatch)

    # compute the labeling on the complete dataset
    labels = km.predict(X)
    assert_equal(v_measure_score(true_labels, labels), 1.0)


def test_k_means_invalid_init():
    k_means = KMeans(init="invalid", n_init=1, k=n_clusters)
    assert_raises(ValueError, k_means.fit, X)


def test_k_means_copyx():
    """Check if copy_x=False returns nearly equal X after de-centering."""
    my_X = X.copy()
    k_means = KMeans(copy_x=False, k=n_clusters, random_state=42).fit(my_X)
    _check_fitted_model(k_means)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)


def test_k_means_singleton():
    """Check k_means with bad initialization and singleton clustering."""
    my_X = np.array([[1.1, 1.1], [0.9, 1.1], [1.1, 0.9], [0.9, 0.9]])
    array_init = np.array([[1.0, 1.0], [5.0, 5.0]])
    k_means = KMeans(init=array_init, k=2, random_state=42, n_init=1)
    k_means.fit(my_X)

    # must be singleton clustering
    assert_equal(np.unique(k_means.labels_).size, 1)


def test_predict():
    k_means = KMeans(k=n_clusters, random_state=42).fit(X)

    # sanity check: predict centroid labels
    pred = k_means.predict(k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = k_means.predict(X)
    assert_array_equal(k_means.predict(X), k_means.labels_)


def test_predict_minibatch_dense_input():
    mb_k_means = MiniBatchKMeans(k=n_clusters, random_state=40).fit(X)

    # sanity check: predict centroid labels
    pred = mb_k_means.predict(mb_k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = mb_k_means.predict(X)
    assert_array_equal(mb_k_means.predict(X), mb_k_means.labels_)


def test_predict_minibatch_sparse_input():
    # TODO: make MiniBatchKMeans more stables using either restarts or
    # multi-init
    mb_k_means = MiniBatchKMeans(k=n_clusters, random_state=40).fit(X_csr)

    # sanity check: re-predict labeling for training set samples
    assert_array_equal(mb_k_means.predict(X_csr), mb_k_means.labels_)

    # sanity check: predict centroid labels
    pred = mb_k_means.predict(mb_k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # check that models trained on sparse input also works for dense input at
    # predict time
    assert_array_equal(mb_k_means.predict(X), mb_k_means.labels_)


def test_input_dtypes():
    X_list = [[0, 0], [10, 10], [12, 9], [-1, 1], [2, 0], [8, 10]]
    X_int = np.array(X_list, dtype=np.int32)
    X_int_csr = sp.csr_matrix(X_int)
    init_int = X_int[:2]

    fitted_models = [
        KMeans(k=2, random_state=42).fit(X_list),
        KMeans(k=2, random_state=42).fit(X_int),
        KMeans(k=2, init=init_int, random_state=42).fit(X_list),
        KMeans(k=2, init=init_int, random_state=42).fit(X_int),
        # mini batch kmeans is very unstable on such a small dataset but we are
        # not interested in testing stability here hence picking up good random
        # init
        MiniBatchKMeans(k=2, random_state=1, chunk_size=2).fit(X_list),
        MiniBatchKMeans(k=2, random_state=1, chunk_size=2).fit(X_int),
        MiniBatchKMeans(k=2, random_state=2, chunk_size=2).fit(X_int_csr),
        MiniBatchKMeans(k=2, random_state=42, init=init_int).fit(X_list),
        MiniBatchKMeans(k=2, random_state=42, init=init_int).fit(X_int),
        MiniBatchKMeans(k=2, random_state=42, init=init_int).fit(X_int_csr),
    ]
    expected_labels = [0, 1, 1, 0, 0, 1]
    scores = np.array([v_measure_score(expected_labels, km.labels_)
                       for km in fitted_models])
    assert_array_equal(scores, np.ones(scores.shape[0]))


def test_transform():
    k_means = KMeans(k=n_clusters)
    k_means.fit(X)
    X_new = k_means.transform(k_means.cluster_centers_)

    for c in range(n_clusters):
        assert_equal(X_new[c, c], 0)
        for c2 in range(n_clusters):
            if c != c2:
                assert_true(X_new[c, c2] > 0)
