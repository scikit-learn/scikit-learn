"""Testing for K-means"""
import sys

import numpy as np
from scipy import sparse as sp

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import if_safe_multiprocessing_with_blas
from sklearn.utils.testing import if_not_mac_os
from sklearn.utils.testing import assert_raise_message


from sklearn.utils.extmath import row_norms
from sklearn.metrics.cluster import v_measure_score
from sklearn.cluster import KMeans, k_means
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster.k_means_ import _labels_inertia
from sklearn.cluster.k_means_ import _mini_batch_step
from sklearn.datasets.samples_generator import make_blobs
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.exceptions import DataConversionWarning
from sklearn.metrics.cluster import homogeneity_score


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


def test_kmeans_dtype():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))
    X = (X * 10).astype(np.uint8)
    km = KMeans(n_init=1).fit(X)
    pred_x = assert_warns(DataConversionWarning, km.predict, X)
    assert_array_equal(km.labels_, pred_x)


def test_elkan_results():
    rnd = np.random.RandomState(0)
    X_normal = rnd.normal(size=(50, 10))
    X_blobs, _ = make_blobs(random_state=0)
    km_full = KMeans(algorithm='full', n_clusters=5, random_state=0, n_init=1)
    km_elkan = KMeans(algorithm='elkan', n_clusters=5,
                      random_state=0, n_init=1)
    for X in [X_normal, X_blobs]:
        km_full.fit(X)
        km_elkan.fit(X)
        assert_array_almost_equal(km_elkan.cluster_centers_, km_full.cluster_centers_)
        assert_array_equal(km_elkan.labels_, km_full.labels_)


def test_labels_assignment_and_inertia():
    # pure numpy implementation as easily auditable reference gold
    # implementation
    rng = np.random.RandomState(42)
    noisy_centers = centers + rng.normal(size=centers.shape)
    labels_gold = - np.ones(n_samples, dtype=np.int)
    mindist = np.empty(n_samples)
    mindist.fill(np.infty)
    for center_id in range(n_clusters):
        dist = np.sum((X - noisy_centers[center_id]) ** 2, axis=1)
        labels_gold[dist < mindist] = center_id
        mindist = np.minimum(dist, mindist)
    inertia_gold = mindist.sum()
    assert_true((mindist >= 0.0).all())
    assert_true((labels_gold != -1).all())

    # perform label assignment using the dense array input
    x_squared_norms = (X ** 2).sum(axis=1)
    labels_array, inertia_array = _labels_inertia(
        X, x_squared_norms, noisy_centers)
    assert_array_almost_equal(inertia_array, inertia_gold)
    assert_array_equal(labels_array, labels_gold)

    # perform label assignment using the sparse CSR input
    x_squared_norms_from_csr = row_norms(X_csr, squared=True)
    labels_csr, inertia_csr = _labels_inertia(
        X_csr, x_squared_norms_from_csr, noisy_centers)
    assert_array_almost_equal(inertia_csr, inertia_gold)
    assert_array_equal(labels_csr, labels_gold)


def test_minibatch_update_consistency():
    # Check that dense and sparse minibatch update give the same results
    rng = np.random.RandomState(42)
    old_centers = centers + rng.normal(size=centers.shape)

    new_centers = old_centers.copy()
    new_centers_csr = old_centers.copy()

    counts = np.zeros(new_centers.shape[0], dtype=np.int32)
    counts_csr = np.zeros(new_centers.shape[0], dtype=np.int32)

    x_squared_norms = (X ** 2).sum(axis=1)
    x_squared_norms_csr = row_norms(X_csr, squared=True)

    buffer = np.zeros(centers.shape[1], dtype=np.double)
    buffer_csr = np.zeros(centers.shape[1], dtype=np.double)

    # extract a small minibatch
    X_mb = X[:10]
    X_mb_csr = X_csr[:10]
    x_mb_squared_norms = x_squared_norms[:10]
    x_mb_squared_norms_csr = x_squared_norms_csr[:10]

    # step 1: compute the dense minibatch update
    old_inertia, incremental_diff = _mini_batch_step(
        X_mb, x_mb_squared_norms, new_centers, counts,
        buffer, 1, None, random_reassign=False)
    assert_greater(old_inertia, 0.0)

    # compute the new inertia on the same batch to check that it decreased
    labels, new_inertia = _labels_inertia(
        X_mb, x_mb_squared_norms, new_centers)
    assert_greater(new_inertia, 0.0)
    assert_less(new_inertia, old_inertia)

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers - old_centers) ** 2)
    assert_almost_equal(incremental_diff, effective_diff)

    # step 2: compute the sparse minibatch update
    old_inertia_csr, incremental_diff_csr = _mini_batch_step(
        X_mb_csr, x_mb_squared_norms_csr, new_centers_csr, counts_csr,
        buffer_csr, 1, None, random_reassign=False)
    assert_greater(old_inertia_csr, 0.0)

    # compute the new inertia on the same batch to check that it decreased
    labels_csr, new_inertia_csr = _labels_inertia(
        X_mb_csr, x_mb_squared_norms_csr, new_centers_csr)
    assert_greater(new_inertia_csr, 0.0)
    assert_less(new_inertia_csr, old_inertia_csr)

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers_csr - old_centers) ** 2)
    assert_almost_equal(incremental_diff_csr, effective_diff)

    # step 3: check that sparse and dense updates lead to the same results
    assert_array_equal(labels, labels_csr)
    assert_array_almost_equal(new_centers, new_centers_csr)
    assert_almost_equal(incremental_diff, incremental_diff_csr)
    assert_almost_equal(old_inertia, old_inertia_csr)
    assert_almost_equal(new_inertia, new_inertia_csr)


def _check_fitted_model(km):
    # check that the number of clusters centers and distinct labels match
    # the expectation
    centers = km.cluster_centers_
    assert_equal(centers.shape, (n_clusters, n_features))

    labels = km.labels_
    assert_equal(np.unique(labels).shape[0], n_clusters)

    # check that the labels assignment are perfect (up to a permutation)
    assert_equal(v_measure_score(true_labels, labels), 1.0)
    assert_greater(km.inertia_, 0.0)

    # check error on dataset being too small
    assert_raises(ValueError, km.fit, [[0., 1.]])


def test_k_means_plus_plus_init():
    km = KMeans(init="k-means++", n_clusters=n_clusters,
                random_state=42).fit(X)
    _check_fitted_model(km)


def test_k_means_new_centers():
    # Explore the part of the code where a new center is reassigned
    X = np.array([[0, 0, 1, 1],
                  [0, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 1, 0, 0]])
    labels = [0, 1, 2, 1, 1, 2]
    bad_centers = np.array([[+0, 1, 0, 0],
                            [.2, 0, .2, .2],
                            [+0, 0, 0, 0]])

    km = KMeans(n_clusters=3, init=bad_centers, n_init=1, max_iter=10,
                random_state=1)
    for this_X in (X, sp.coo_matrix(X)):
        km.fit(this_X)
        this_labels = km.labels_
        # Reorder the labels so that the first instance is in cluster 0,
        # the second in cluster 1, ...
        this_labels = np.unique(this_labels, return_index=True)[1][this_labels]
        np.testing.assert_array_equal(this_labels, labels)


@if_safe_multiprocessing_with_blas
def test_k_means_plus_plus_init_2_jobs():
    if sys.version_info[:2] < (3, 4):
        raise SkipTest(
            "Possible multi-process bug with some BLAS under Python < 3.4")

    km = KMeans(init="k-means++", n_clusters=n_clusters, n_jobs=2,
                random_state=42).fit(X)
    _check_fitted_model(km)


def test_k_means_precompute_distances_flag():
    # check that a warning is raised if the precompute_distances flag is not
    # supported
    km = KMeans(precompute_distances="wrong")
    assert_raises(ValueError, km.fit, X)


def test_k_means_plus_plus_init_sparse():
    km = KMeans(init="k-means++", n_clusters=n_clusters, random_state=42)
    km.fit(X_csr)
    _check_fitted_model(km)


def test_k_means_random_init():
    km = KMeans(init="random", n_clusters=n_clusters, random_state=42)
    km.fit(X)
    _check_fitted_model(km)


def test_k_means_random_init_sparse():
    km = KMeans(init="random", n_clusters=n_clusters, random_state=42)
    km.fit(X_csr)
    _check_fitted_model(km)


def test_k_means_plus_plus_init_not_precomputed():
    km = KMeans(init="k-means++", n_clusters=n_clusters, random_state=42,
                precompute_distances=False).fit(X)
    _check_fitted_model(km)


def test_k_means_random_init_not_precomputed():
    km = KMeans(init="random", n_clusters=n_clusters, random_state=42,
                precompute_distances=False).fit(X)
    _check_fitted_model(km)


def test_k_means_perfect_init():
    km = KMeans(init=centers.copy(), n_clusters=n_clusters, random_state=42,
                n_init=1)
    km.fit(X)
    _check_fitted_model(km)


def test_k_means_n_init():
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 2))

    # two regression tests on bad n_init argument
    # previous bug: n_init <= 0 threw non-informative TypeError (#3858)
    assert_raises_regex(ValueError, "n_init", KMeans(n_init=0).fit, X)
    assert_raises_regex(ValueError, "n_init", KMeans(n_init=-1).fit, X)


def test_k_means_explicit_init_shape():
    # test for sensible errors when giving explicit init
    # with wrong number of features or clusters
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(40, 3))
    for Class in [KMeans, MiniBatchKMeans]:
        # mismatch of number of features
        km = Class(n_init=1, init=X[:, :2], n_clusters=len(X))
        msg = "does not match the number of features of the data"
        assert_raises_regex(ValueError, msg, km.fit, X)
        # for callable init
        km = Class(n_init=1, init=lambda X_, k, random_state: X_[:, :2], n_clusters=len(X))
        assert_raises_regex(ValueError, msg, km.fit, X)
        # mismatch of number of clusters
        msg = "does not match the number of clusters"
        km = Class(n_init=1, init=X[:2, :], n_clusters=3)
        assert_raises_regex(ValueError, msg, km.fit, X)
        # for callable init
        km = Class(n_init=1, init=lambda X_, k, random_state: X_[:2, :], n_clusters=3)
        assert_raises_regex(ValueError, msg, km.fit, X)


def test_k_means_fortran_aligned_data():
    # Check the KMeans will work well, even if X is a fortran-aligned data.
    X = np.asfortranarray([[0, 0], [0, 1], [0, 1]])
    centers = np.array([[0, 0], [0, 1]])
    labels = np.array([0, 1, 1])
    km = KMeans(n_init=1, init=centers, precompute_distances=False,
                random_state=42, n_clusters=2)
    km.fit(X)
    assert_array_equal(km.cluster_centers_, centers)
    assert_array_equal(km.labels_, labels)


def test_mb_k_means_plus_plus_init_dense_array():
    mb_k_means = MiniBatchKMeans(init="k-means++", n_clusters=n_clusters,
                                 random_state=42)
    mb_k_means.fit(X)
    _check_fitted_model(mb_k_means)


def test_mb_kmeans_verbose():
    mb_k_means = MiniBatchKMeans(init="k-means++", n_clusters=n_clusters,
                                 random_state=42, verbose=1)
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        mb_k_means.fit(X)
    finally:
        sys.stdout = old_stdout


def test_mb_k_means_plus_plus_init_sparse_matrix():
    mb_k_means = MiniBatchKMeans(init="k-means++", n_clusters=n_clusters,
                                 random_state=42)
    mb_k_means.fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_minibatch_init_with_large_k():
    mb_k_means = MiniBatchKMeans(init='k-means++', init_size=10, n_clusters=20)
    # Check that a warning is raised, as the number clusters is larger
    # than the init_size
    assert_warns(RuntimeWarning, mb_k_means.fit, X)


def test_minibatch_k_means_random_init_dense_array():
    # increase n_init to make random init stable enough
    mb_k_means = MiniBatchKMeans(init="random", n_clusters=n_clusters,
                                 random_state=42, n_init=10).fit(X)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_random_init_sparse_csr():
    # increase n_init to make random init stable enough
    mb_k_means = MiniBatchKMeans(init="random", n_clusters=n_clusters,
                                 random_state=42, n_init=10).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_perfect_init_dense_array():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), n_clusters=n_clusters,
                                 random_state=42, n_init=1).fit(X)
    _check_fitted_model(mb_k_means)


def test_minibatch_k_means_init_multiple_runs_with_explicit_centers():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), n_clusters=n_clusters,
                                 random_state=42, n_init=10)
    assert_warns(RuntimeWarning, mb_k_means.fit, X)


def test_minibatch_k_means_perfect_init_sparse_csr():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), n_clusters=n_clusters,
                                 random_state=42, n_init=1).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_minibatch_sensible_reassign_fit():
    # check if identical initial clusters are reassigned
    # also a regression test for when there are more desired reassignments than
    # samples.
    zeroed_X, true_labels = make_blobs(n_samples=100, centers=5,
                                       cluster_std=1., random_state=42)
    zeroed_X[::2, :] = 0
    mb_k_means = MiniBatchKMeans(n_clusters=20, batch_size=10, random_state=42,
                                 init="random")
    mb_k_means.fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert_greater(mb_k_means.cluster_centers_.any(axis=1).sum(), 10)

    # do the same with batch-size > X.shape[0] (regression test)
    mb_k_means = MiniBatchKMeans(n_clusters=20, batch_size=201,
                                 random_state=42, init="random")
    mb_k_means.fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert_greater(mb_k_means.cluster_centers_.any(axis=1).sum(), 10)


def test_minibatch_sensible_reassign_partial_fit():
    zeroed_X, true_labels = make_blobs(n_samples=n_samples, centers=5,
                                       cluster_std=1., random_state=42)
    zeroed_X[::2, :] = 0
    mb_k_means = MiniBatchKMeans(n_clusters=20, random_state=42, init="random")
    for i in range(100):
        mb_k_means.partial_fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert_greater(mb_k_means.cluster_centers_.any(axis=1).sum(), 10)


def test_minibatch_reassign():
    # Give a perfect initialization, but a large reassignment_ratio,
    # as a result all the centers should be reassigned and the model
    # should not longer be good
    for this_X in (X, X_csr):
        mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, batch_size=100,
                                     random_state=42)
        mb_k_means.fit(this_X)

        score_before = mb_k_means.score(this_X)
        try:
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            # Turn on verbosity to smoke test the display code
            _mini_batch_step(this_X, (X ** 2).sum(axis=1),
                             mb_k_means.cluster_centers_,
                             mb_k_means.counts_,
                             np.zeros(X.shape[1], np.double),
                             False, distances=np.zeros(X.shape[0]),
                             random_reassign=True, random_state=42,
                             reassignment_ratio=1, verbose=True)
        finally:
            sys.stdout = old_stdout
        assert_greater(score_before, mb_k_means.score(this_X))

    # Give a perfect initialization, with a small reassignment_ratio,
    # no center should be reassigned
    for this_X in (X, X_csr):
        mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, batch_size=100,
                                     init=centers.copy(),
                                     random_state=42, n_init=1)
        mb_k_means.fit(this_X)
        clusters_before = mb_k_means.cluster_centers_
        # Turn on verbosity to smoke test the display code
        _mini_batch_step(this_X, (X ** 2).sum(axis=1),
                         mb_k_means.cluster_centers_,
                         mb_k_means.counts_,
                         np.zeros(X.shape[1], np.double),
                         False, distances=np.zeros(X.shape[0]),
                         random_reassign=True, random_state=42,
                         reassignment_ratio=1e-15)
        assert_array_almost_equal(clusters_before, mb_k_means.cluster_centers_)


def test_minibatch_with_many_reassignments():
    # Test for the case that the number of clusters to reassign is bigger
    # than the batch_size
    n_samples = 550
    rnd = np.random.RandomState(42)
    X = rnd.uniform(size=(n_samples, 10))
    # Check that the fit works if n_clusters is bigger than the batch_size.
    # Run the test with 550 clusters and 550 samples, because it turned out
    # that this values ensure that the number of clusters to reassign
    # is always bigger than the batch_size
    n_clusters = 550
    MiniBatchKMeans(n_clusters=n_clusters,
                    batch_size=100,
                    init_size=n_samples,
                    random_state=42).fit(X)


def test_sparse_mb_k_means_callable_init():

    def test_init(X, k, random_state):
        return centers

    # Small test to check that giving the wrong number of centers
    # raises a meaningful error
    msg = "does not match the number of clusters"
    assert_raises_regex(ValueError, msg, MiniBatchKMeans(init=test_init,
                                                         random_state=42).fit,
                        X_csr)

    # Now check that the fit actually works
    mb_k_means = MiniBatchKMeans(n_clusters=3, init=test_init,
                                 random_state=42).fit(X_csr)
    _check_fitted_model(mb_k_means)


def test_mini_batch_k_means_random_init_partial_fit():
    km = MiniBatchKMeans(n_clusters=n_clusters, init="random", random_state=42)

    # use the partial_fit API for online learning
    for X_minibatch in np.array_split(X, 10):
        km.partial_fit(X_minibatch)

    # compute the labeling on the complete dataset
    labels = km.predict(X)
    assert_equal(v_measure_score(true_labels, labels), 1.0)


def test_minibatch_default_init_size():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), n_clusters=n_clusters,
                                 batch_size=10, random_state=42,
                                 n_init=1).fit(X)
    assert_equal(mb_k_means.init_size_, 3 * mb_k_means.batch_size)
    _check_fitted_model(mb_k_means)


def test_minibatch_tol():
    mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, batch_size=10,
                                 random_state=42, tol=.01).fit(X)
    _check_fitted_model(mb_k_means)


def test_minibatch_set_init_size():
    mb_k_means = MiniBatchKMeans(init=centers.copy(), n_clusters=n_clusters,
                                 init_size=666, random_state=42,
                                 n_init=1).fit(X)
    assert_equal(mb_k_means.init_size, 666)
    assert_equal(mb_k_means.init_size_, n_samples)
    _check_fitted_model(mb_k_means)


def test_k_means_invalid_init():
    km = KMeans(init="invalid", n_init=1, n_clusters=n_clusters)
    assert_raises(ValueError, km.fit, X)


def test_mini_match_k_means_invalid_init():
    km = MiniBatchKMeans(init="invalid", n_init=1, n_clusters=n_clusters)
    assert_raises(ValueError, km.fit, X)


def test_k_means_copyx():
    # Check if copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    km = KMeans(copy_x=False, n_clusters=n_clusters, random_state=42)
    km.fit(my_X)
    _check_fitted_model(km)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)


def test_k_means_non_collapsed():
    # Check k_means with a bad initialization does not yield a singleton
    # Starting with bad centers that are quickly ignored should not
    # result in a repositioning of the centers to the center of mass that
    # would lead to collapsed centers which in turns make the clustering
    # dependent of the numerical unstabilities.
    my_X = np.array([[1.1, 1.1], [0.9, 1.1], [1.1, 0.9], [0.9, 1.1]])
    array_init = np.array([[1.0, 1.0], [5.0, 5.0], [-5.0, -5.0]])
    km = KMeans(init=array_init, n_clusters=3, random_state=42, n_init=1)
    km.fit(my_X)

    # centers must not been collapsed
    assert_equal(len(np.unique(km.labels_)), 3)

    centers = km.cluster_centers_
    assert_true(np.linalg.norm(centers[0] - centers[1]) >= 0.1)
    assert_true(np.linalg.norm(centers[0] - centers[2]) >= 0.1)
    assert_true(np.linalg.norm(centers[1] - centers[2]) >= 0.1)


def test_predict():
    km = KMeans(n_clusters=n_clusters, random_state=42)

    km.fit(X)

    # sanity check: predict centroid labels
    pred = km.predict(km.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = km.predict(X)
    assert_array_equal(pred, km.labels_)

    # re-predict labels for training set using fit_predict
    pred = km.fit_predict(X)
    assert_array_equal(pred, km.labels_)


def test_score():

    km1 = KMeans(n_clusters=n_clusters, max_iter=1, random_state=42, n_init=1)
    s1 = km1.fit(X).score(X)
    km2 = KMeans(n_clusters=n_clusters, max_iter=10, random_state=42, n_init=1)
    s2 = km2.fit(X).score(X)
    assert_greater(s2, s1)

    km1 = KMeans(n_clusters=n_clusters, max_iter=1, random_state=42, n_init=1,
                 algorithm='elkan')
    s1 = km1.fit(X).score(X)
    km2 = KMeans(n_clusters=n_clusters, max_iter=10, random_state=42, n_init=1,
                 algorithm='elkan')
    s2 = km2.fit(X).score(X)
    assert_greater(s2, s1)


def test_predict_minibatch_dense_input():
    mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, random_state=40).fit(X)

    # sanity check: predict centroid labels
    pred = mb_k_means.predict(mb_k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # sanity check: re-predict labeling for training set samples
    pred = mb_k_means.predict(X)
    assert_array_equal(mb_k_means.predict(X), mb_k_means.labels_)


def test_predict_minibatch_kmeanspp_init_sparse_input():
    mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, init='k-means++',
                                 n_init=10).fit(X_csr)

    # sanity check: re-predict labeling for training set samples
    assert_array_equal(mb_k_means.predict(X_csr), mb_k_means.labels_)

    # sanity check: predict centroid labels
    pred = mb_k_means.predict(mb_k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # check that models trained on sparse input also works for dense input at
    # predict time
    assert_array_equal(mb_k_means.predict(X), mb_k_means.labels_)


def test_predict_minibatch_random_init_sparse_input():
    mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, init='random',
                                 n_init=10).fit(X_csr)

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
        KMeans(n_clusters=2).fit(X_list),
        KMeans(n_clusters=2).fit(X_int),
        KMeans(n_clusters=2, init=init_int, n_init=1).fit(X_list),
        KMeans(n_clusters=2, init=init_int, n_init=1).fit(X_int),
        # mini batch kmeans is very unstable on such a small dataset hence
        # we use many inits
        MiniBatchKMeans(n_clusters=2, n_init=10, batch_size=2).fit(X_list),
        MiniBatchKMeans(n_clusters=2, n_init=10, batch_size=2).fit(X_int),
        MiniBatchKMeans(n_clusters=2, n_init=10, batch_size=2).fit(X_int_csr),
        MiniBatchKMeans(n_clusters=2, batch_size=2,
                        init=init_int, n_init=1).fit(X_list),
        MiniBatchKMeans(n_clusters=2, batch_size=2,
                        init=init_int, n_init=1).fit(X_int),
        MiniBatchKMeans(n_clusters=2, batch_size=2,
                        init=init_int, n_init=1).fit(X_int_csr),
    ]
    expected_labels = [0, 1, 1, 0, 0, 1]
    scores = np.array([v_measure_score(expected_labels, km.labels_)
                       for km in fitted_models])
    assert_array_equal(scores, np.ones(scores.shape[0]))


def test_transform():
    km = KMeans(n_clusters=n_clusters)
    km.fit(X)
    X_new = km.transform(km.cluster_centers_)

    for c in range(n_clusters):
        assert_equal(X_new[c, c], 0)
        for c2 in range(n_clusters):
            if c != c2:
                assert_greater(X_new[c, c2], 0)


def test_fit_transform():
    X1 = KMeans(n_clusters=3, random_state=51).fit(X).transform(X)
    X2 = KMeans(n_clusters=3, random_state=51).fit_transform(X)
    assert_array_equal(X1, X2)


def test_predict_equal_labels():
    km = KMeans(random_state=13, n_jobs=1, n_init=1, max_iter=1,
                algorithm='full')
    km.fit(X)
    assert_array_equal(km.predict(X), km.labels_)

    km = KMeans(random_state=13, n_jobs=1, n_init=1, max_iter=1,
                algorithm='elkan')
    km.fit(X)
    assert_array_equal(km.predict(X), km.labels_)


def test_full_vs_elkan():

    km1 = KMeans(algorithm='full', random_state=13)
    km2 = KMeans(algorithm='elkan', random_state=13)

    km1.fit(X)
    km2.fit(X)

    homogeneity_score(km1.predict(X), km2.predict(X)) == 1.0


def test_n_init():
    # Check that increasing the number of init increases the quality
    n_runs = 5
    n_init_range = [1, 5, 10]
    inertia = np.zeros((len(n_init_range), n_runs))
    for i, n_init in enumerate(n_init_range):
        for j in range(n_runs):
            km = KMeans(n_clusters=n_clusters, init="random", n_init=n_init,
                        random_state=j).fit(X)
            inertia[i, j] = km.inertia_

    inertia = inertia.mean(axis=1)
    failure_msg = ("Inertia %r should be decreasing"
                   " when n_init is increasing.") % list(inertia)
    for i in range(len(n_init_range) - 1):
        assert_true(inertia[i] >= inertia[i + 1], failure_msg)


def test_k_means_function():
    # test calling the k_means function directly
    # catch output
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        cluster_centers, labels, inertia = k_means(X, n_clusters=n_clusters,
                                                   verbose=True)
    finally:
        sys.stdout = old_stdout
    centers = cluster_centers
    assert_equal(centers.shape, (n_clusters, n_features))

    labels = labels
    assert_equal(np.unique(labels).shape[0], n_clusters)

    # check that the labels assignment are perfect (up to a permutation)
    assert_equal(v_measure_score(true_labels, labels), 1.0)
    assert_greater(inertia, 0.0)

    # check warning when centers are passed
    assert_warns(RuntimeWarning, k_means, X, n_clusters=n_clusters,
                 init=centers)

    # to many clusters desired
    assert_raises(ValueError, k_means, X, n_clusters=X.shape[0] + 1)


def test_x_squared_norms_init_centroids():
    """Test that x_squared_norms can be None in _init_centroids"""
    from sklearn.cluster.k_means_ import _init_centroids

    X_norms = np.sum(X**2, axis=1)
    precompute = _init_centroids(
        X, 3, "k-means++", random_state=0, x_squared_norms=X_norms)
    assert_array_equal(
        precompute,
        _init_centroids(X, 3, "k-means++", random_state=0))


def test_max_iter_error():

    km = KMeans(max_iter=-1)
    assert_raise_message(ValueError, 'Number of iterations should be', km.fit, X)
