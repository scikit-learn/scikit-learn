"""Testing for K-means"""
import sys

import numpy as np
from scipy import sparse as sp
from threadpoolctl import threadpool_limits

import pytest

from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import assert_almost_equal
from sklearn.utils.fixes import _astype_copy_false
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning

from sklearn.utils.extmath import row_norms
from sklearn.metrics import pairwise_distances
from sklearn.metrics import pairwise_distances_argmin
from sklearn.metrics.cluster import v_measure_score
from sklearn.cluster import KMeans, k_means
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster._kmeans import _labels_inertia
from sklearn.cluster._kmeans import _mini_batch_step
from sklearn.cluster._k_means_fast import _relocate_empty_clusters_dense
from sklearn.cluster._k_means_fast import _relocate_empty_clusters_sparse
from sklearn.cluster._k_means_fast import _euclidean_dense_dense_wrapper
from sklearn.cluster._k_means_fast import _euclidean_sparse_dense_wrapper
from sklearn.cluster._k_means_fast import _inertia_dense
from sklearn.cluster._k_means_fast import _inertia_sparse
from sklearn.datasets import make_blobs
from io import StringIO


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


@pytest.mark.parametrize("array_constr", [np.array, sp.csr_matrix],
                         ids=["dense", "sparse"])
@pytest.mark.parametrize("algo", ["full", "elkan"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_kmeans_results(array_constr, algo, dtype):
    # Checks that KMeans works as intended on toy dataset by comparing with
    # expected results computed by hand.
    X = array_constr([[0, 0], [0.5, 0], [0.5, 1], [1, 1]], dtype=dtype)
    sample_weight = [3, 1, 1, 3]
    init_centers = np.array([[0, 0], [1, 1]], dtype=dtype)

    expected_labels = [0, 0, 1, 1]
    expected_inertia = 0.375
    expected_centers = np.array([[0.125, 0], [0.875, 1]], dtype=dtype)
    expected_n_iter = 2

    kmeans = KMeans(n_clusters=2, n_init=1, init=init_centers, algorithm=algo)
    kmeans.fit(X, sample_weight=sample_weight)

    assert_array_equal(kmeans.labels_, expected_labels)
    assert_allclose(kmeans.inertia_, expected_inertia)
    assert_allclose(kmeans.cluster_centers_, expected_centers)
    assert kmeans.n_iter_ == expected_n_iter


@pytest.mark.parametrize("array_constr", [np.array, sp.csr_matrix],
                         ids=['dense', 'sparse'])
@pytest.mark.parametrize("algo", ['full', 'elkan'])
def test_kmeans_relocated_clusters(array_constr, algo):
    # check that empty clusters are relocated as expected
    X = array_constr([[0, 0], [0.5, 0], [0.5, 1], [1, 1]])

    # second center too far from others points will be empty at first iter
    init_centers = np.array([[0.5, 0.5], [3, 3]])

    expected_labels = [0, 0, 1, 1]
    expected_inertia = 0.25
    expected_centers = [[0.25, 0], [0.75, 1]]
    expected_n_iter = 3

    kmeans = KMeans(n_clusters=2, n_init=1, init=init_centers, algorithm=algo)
    kmeans.fit(X)

    assert_array_equal(kmeans.labels_, expected_labels)
    assert_allclose(kmeans.inertia_, expected_inertia)
    assert_allclose(kmeans.cluster_centers_, expected_centers)
    assert kmeans.n_iter_ == expected_n_iter


@pytest.mark.parametrize("array_constr", [np.array, sp.csr_matrix],
                         ids=["dense", "sparse"])
def test_relocate_empty_clusters(array_constr):
    # test for the _relocate_empty_clusters_(dense/sparse) helpers

    # Synthetic dataset with 3 obvious clusters of different sizes
    X = np.array(
        [-10., -9.5, -9, -8.5, -8, -1, 1, 9, 9.5, 10]).reshape(-1, 1)
    X = array_constr(X)
    sample_weight = np.ones(10)

    # centers all initialized to the first point of X
    centers_old = np.array([-10., -10, -10]).reshape(-1, 1)

    # With this initialization, all points will be assigned to the first center
    # At this point a center in centers_new is the weighted sum of the points
    # it contains if it's not empty, otherwise it is the same as before.
    centers_new = np.array([-16.5, -10, -10]).reshape(-1, 1)
    weight_in_clusters = np.array([10., 0, 0])
    labels = np.zeros(10, dtype=np.int32)

    if array_constr is np.array:
        _relocate_empty_clusters_dense(X, sample_weight, centers_old,
                                       centers_new, weight_in_clusters, labels)
    else:
        _relocate_empty_clusters_sparse(X.data, X.indices, X.indptr,
                                        sample_weight, centers_old,
                                        centers_new, weight_in_clusters,
                                        labels)

    # The relocation scheme will take the 2 points farthest from the center and
    # assign them to the 2 empty clusters, i.e. points at 10 and at 9.9. The
    # first center will be updated to contain the other 8 points.
    assert_array_equal(weight_in_clusters, [8, 1, 1])
    assert_allclose(centers_new, [[-36], [10], [9.5]])


@pytest.mark.parametrize("distribution", ["normal", "blobs"])
@pytest.mark.parametrize("array_constr", [np.array, sp.csr_matrix],
                         ids=["dense", "sparse"])
@pytest.mark.parametrize("tol", [1e-2, 1e-4, 1e-8])
def test_kmeans_elkan_results(distribution, array_constr, tol):
    # Check that results are identical between lloyd and elkan algorithms
    rnd = np.random.RandomState(0)
    if distribution == "normal":
        X = rnd.normal(size=(5000, 10))
    else:
        X, _ = make_blobs(random_state=rnd)
    X[X < 0] = 0
    X = array_constr(X)

    km_full = KMeans(algorithm="full", n_clusters=5,
                     random_state=0, n_init=1, tol=tol)
    km_elkan = KMeans(algorithm="elkan", n_clusters=5,
                      random_state=0, n_init=1, tol=tol)

    km_full.fit(X)
    km_elkan.fit(X)
    assert_allclose(km_elkan.cluster_centers_, km_full.cluster_centers_)
    assert_array_equal(km_elkan.labels_, km_full.labels_)
    assert km_elkan.n_iter_ == km_full.n_iter_
    assert km_elkan.inertia_ == pytest.approx(km_full.inertia_, rel=1e-6)


@pytest.mark.parametrize("algorithm", ["full", "elkan"])
def test_kmeans_convergence(algorithm):
    # Check that KMeans stops when convergence is reached when tol=0. (#16075)
    # We can only ensure that if the number of threads is not to large,
    # otherwise the roundings errors coming from the unpredictability of
    # the order in which chunks are processed make the convergence criterion
    # to never be exactly 0.
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(5000, 10))

    with threadpool_limits(limits=1, user_api="openmp"):
        km = KMeans(algorithm=algorithm, n_clusters=5, random_state=0,
                    n_init=1, tol=0, max_iter=300).fit(X)

    assert km.n_iter_ < 300


def test_minibatch_update_consistency():
    # Check that dense and sparse minibatch update give the same results
    rng = np.random.RandomState(42)
    old_centers = centers + rng.normal(size=centers.shape)

    new_centers = old_centers.copy()
    new_centers_csr = old_centers.copy()

    weight_sums = np.zeros(new_centers.shape[0], dtype=np.double)
    weight_sums_csr = np.zeros(new_centers.shape[0], dtype=np.double)

    x_squared_norms = (X ** 2).sum(axis=1)
    x_squared_norms_csr = row_norms(X_csr, squared=True)

    buffer = np.zeros(centers.shape[1], dtype=np.double)
    buffer_csr = np.zeros(centers.shape[1], dtype=np.double)

    # extract a small minibatch
    X_mb = X[:10]
    X_mb_csr = X_csr[:10]
    x_mb_squared_norms = x_squared_norms[:10]
    x_mb_squared_norms_csr = x_squared_norms_csr[:10]

    sample_weight_mb = np.ones(X_mb.shape[0], dtype=np.double)

    # step 1: compute the dense minibatch update
    old_inertia, incremental_diff = _mini_batch_step(
        X_mb, sample_weight_mb, x_mb_squared_norms, new_centers, weight_sums,
        buffer, 1, None, random_reassign=False)
    assert old_inertia > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels, new_inertia = _labels_inertia(
        X_mb, sample_weight_mb, x_mb_squared_norms, new_centers)
    assert new_inertia > 0.0
    assert new_inertia < old_inertia

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers - old_centers) ** 2)
    assert_almost_equal(incremental_diff, effective_diff)

    # step 2: compute the sparse minibatch update
    old_inertia_csr, incremental_diff_csr = _mini_batch_step(
        X_mb_csr, sample_weight_mb, x_mb_squared_norms_csr, new_centers_csr,
        weight_sums_csr, buffer_csr, 1, None, random_reassign=False)
    assert old_inertia_csr > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels_csr, new_inertia_csr = _labels_inertia(
        X_mb_csr, sample_weight_mb, x_mb_squared_norms_csr, new_centers_csr)
    assert new_inertia_csr > 0.0
    assert new_inertia_csr < old_inertia_csr

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
    assert centers.shape == (n_clusters, n_features)

    labels = km.labels_
    assert np.unique(labels).shape[0] == n_clusters

    # check that the labels assignment are perfect (up to a permutation)
    assert v_measure_score(true_labels, labels) == 1.0
    assert km.inertia_ > 0.0


@pytest.mark.parametrize("data", [X, X_csr], ids=["dense", "sparse"])
@pytest.mark.parametrize("init", ["random", "k-means++", centers,
                                  lambda X, k, random_state: centers],
                         ids=["random", "k-means++", "ndarray", "callable"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_all_init(Estimator, data, init):
    # Check KMeans and MiniBatchKMeans with all possible init.
    n_init = 10 if type(init) is str else 1
    km = Estimator(init=init, n_clusters=n_clusters, random_state=42,
                   n_init=n_init).fit(data)
    _check_fitted_model(km)


@pytest.mark.parametrize("init", ["random", "k-means++", centers,
                                  lambda X, k, random_state: centers],
                         ids=["random", "k-means++", "ndarray", "callable"])
def test_minibatch_kmeans_partial_fit_init(init):
    # Check MiniBatchKMeans init with partial_fit
    km = MiniBatchKMeans(init=init, n_clusters=n_clusters, random_state=0)
    for i in range(100):
        # "random" init requires many batches to recover the true labels.
        km.partial_fit(X)
    _check_fitted_model(km)


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_fortran_aligned_data(Estimator):
    # Check that KMeans works with fortran-aligned data.
    X_fortran = np.asfortranarray(X)
    centers_fortran = np.asfortranarray(centers)

    km_c = Estimator(n_clusters=n_clusters, init=centers, n_init=1,
                     random_state=42).fit(X)
    km_f = Estimator(n_clusters=n_clusters, init=centers_fortran, n_init=1,
                     random_state=42).fit(X_fortran)
    assert_allclose(km_c.cluster_centers_, km_f.cluster_centers_)
    assert_array_equal(km_c.labels_, km_f.labels_)


@pytest.mark.parametrize('algo', ['full', 'elkan'])
@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('constructor', [np.asarray, sp.csr_matrix])
@pytest.mark.parametrize('seed, max_iter, tol', [
    (0, 2, 1e-7),    # strict non-convergence
    (1, 2, 1e-1),    # loose non-convergence
    (3, 300, 1e-7),  # strict convergence
    (4, 300, 1e-1),  # loose convergence
])
def test_k_means_fit_predict(algo, dtype, constructor, seed, max_iter, tol):
    # check that fit.predict gives same result as fit_predict
    # There's a very small chance of failure with elkan on unstructured dataset
    # because predict method uses fast euclidean distances computation which
    # may cause small numerical instabilities.
    # NB: This test is largely redundant with respect to test_predict and
    #     test_predict_equal_labels.  This test has the added effect of
    #     testing idempotence of the fittng procesdure which appears to
    #     be where it fails on some MacOS setups.
    if sys.platform == "darwin":
        pytest.xfail(
            "Known failures on MacOS, See "
            "https://github.com/scikit-learn/scikit-learn/issues/12644")

    rng = np.random.RandomState(seed)

    X = make_blobs(n_samples=1000, n_features=10, centers=10,
                   random_state=rng)[0].astype(dtype, copy=False)
    X = constructor(X)

    kmeans = KMeans(algorithm=algo, n_clusters=10, random_state=seed,
                    tol=tol, max_iter=max_iter)

    labels_1 = kmeans.fit(X).predict(X)
    labels_2 = kmeans.fit_predict(X)

    # Due to randomness in the order in which chunks of data are processed when
    # using more than one thread, the absolute values of the labels can be
    # different between the 2 strategies but they should correspond to the same
    # clustering.
    assert v_measure_score(labels_1, labels_2) == 1


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_verbose(Estimator):
    # Check verbose mode of KMeans and MiniBatchKMeans for better coverage.
    km = Estimator(n_clusters=n_clusters, random_state=42, verbose=1)
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        km.fit(X)
    finally:
        sys.stdout = old_stdout


def test_minibatch_kmeans_warning_init_size():
    # Check that a warning is raised when init_size is smaller than n_clusters
    with pytest.warns(RuntimeWarning,
                      match=r"init_size.* should be larger than n_clusters"):
        MiniBatchKMeans(init_size=10, n_clusters=20).fit(X)


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_warning_n_init_precomputed_centers(Estimator):
    # Check that a warning is raised when n_init > 1 and an array is passed for
    # the init parameter.
    with pytest.warns(RuntimeWarning,
                      match="Explicit initial center position passed: "
                            "performing only one init"):
        Estimator(init=centers, n_clusters=n_clusters, n_init=10).fit(X)


def test_minibatch_sensible_reassign():
    # check that identical initial clusters are reassigned
    # also a regression test for when there are more desired reassignments than
    # samples.
    zeroed_X, true_labels = make_blobs(n_samples=100, centers=5,
                                       random_state=42)
    zeroed_X[::2, :] = 0

    km = MiniBatchKMeans(n_clusters=20, batch_size=10, random_state=42,
                         init="random").fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert km.cluster_centers_.any(axis=1).sum() > 10

    # do the same with batch-size > X.shape[0] (regression test)
    km = MiniBatchKMeans(n_clusters=20, batch_size=200, random_state=42,
                         init="random").fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert km.cluster_centers_.any(axis=1).sum() > 10

    # do the same with partial_fit API
    km = MiniBatchKMeans(n_clusters=20, random_state=42, init="random")
    for i in range(100):
        km.partial_fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    assert km.cluster_centers_.any(axis=1).sum() > 10


def test_minibatch_reassign():
    # Give a perfect initialization, but a large reassignment_ratio,
    # as a result all the centers should be reassigned and the model
    # should no longer be good
    sample_weight = np.ones(X.shape[0], dtype=X.dtype)
    for this_X in (X, X_csr):
        mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, batch_size=100,
                                     random_state=42)
        mb_k_means.fit(this_X)

        score_before = mb_k_means.score(this_X)
        try:
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            # Turn on verbosity to smoke test the display code
            _mini_batch_step(this_X, sample_weight, (X ** 2).sum(axis=1),
                             mb_k_means.cluster_centers_,
                             mb_k_means._counts,
                             np.zeros(X.shape[1], np.double),
                             False, distances=np.zeros(X.shape[0]),
                             random_reassign=True, random_state=42,
                             reassignment_ratio=1, verbose=True)
        finally:
            sys.stdout = old_stdout
        assert score_before > mb_k_means.score(this_X)

    # Give a perfect initialization, with a small reassignment_ratio,
    # no center should be reassigned
    for this_X in (X, X_csr):
        mb_k_means = MiniBatchKMeans(n_clusters=n_clusters, batch_size=100,
                                     init=centers.copy(),
                                     random_state=42, n_init=1)
        mb_k_means.fit(this_X)
        clusters_before = mb_k_means.cluster_centers_
        # Turn on verbosity to smoke test the display code
        _mini_batch_step(this_X, sample_weight, (X ** 2).sum(axis=1),
                         mb_k_means.cluster_centers_,
                         mb_k_means._counts,
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
    assert v_measure_score(true_labels, labels) == 1.0


def test_minibatch_kmeans_init_size():
    # Check the internal _init_size attribute of MiniBatchKMeans

    # default init size should be 3 * batch_size
    km = MiniBatchKMeans(n_clusters=10, batch_size=5, n_init=1).fit(X)
    assert km._init_size == 15

    # if 3 * batch size < n_clusters, it should then be 3 * n_clusters
    km = MiniBatchKMeans(n_clusters=10, batch_size=1, n_init=1).fit(X)
    assert km._init_size == 30

    # it should not be larger than n_samples
    km = MiniBatchKMeans(n_clusters=10, batch_size=5, n_init=1,
                         init_size=n_samples + 1).fit(X)
    assert km._init_size == n_samples


def test_kmeans_copyx():
    # Check that copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    km = KMeans(copy_x=False, n_clusters=n_clusters, random_state=42)
    km.fit(my_X)
    _check_fitted_model(km)

    # check that my_X is de-centered
    assert_allclose(my_X, X)


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_score_max_iter(Estimator):
    # Check that fitting KMeans or MiniBatchKMeans with more iterations gives
    # better score
    X = np.random.RandomState(0).randn(100, 10)

    km1 = Estimator(n_init=1, random_state=42, max_iter=1)
    s1 = km1.fit(X).score(X)
    km2 = Estimator(n_init=1, random_state=42, max_iter=10)
    s2 = km2.fit(X).score(X)
    assert s2 > s1


@pytest.mark.parametrize('Estimator', [KMeans, MiniBatchKMeans])
@pytest.mark.parametrize('data', [X, X_csr], ids=['dense', 'sparse'])
@pytest.mark.parametrize('init', ['random', 'k-means++', centers.copy()])
def test_predict(Estimator, data, init):
    n_init = 10 if type(init) is str else 1
    k_means = Estimator(n_clusters=n_clusters, init=init,
                        n_init=n_init, random_state=0).fit(data)

    # sanity check: re-predict labeling for training set samples
    assert_array_equal(k_means.predict(data), k_means.labels_)

    # sanity check: predict centroid labels
    pred = k_means.predict(k_means.cluster_centers_)
    assert_array_equal(pred, np.arange(n_clusters))

    # re-predict labels for training set using fit_predict
    pred = k_means.fit_predict(data)
    assert_array_equal(pred, k_means.labels_)


@pytest.mark.parametrize("init", ["random", "k-means++", centers],
                         ids=["random", "k-means++", "ndarray"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_predict_dense_sparse(Estimator, init):
    # check that models trained on sparse input also works for dense input at
    # predict time and vice versa.
    n_init = 10 if type(init) is str else 1
    km = Estimator(n_clusters=n_clusters, init=init, n_init=n_init,
                   random_state=0)

    km.fit(X_csr)
    assert_array_equal(km.predict(X), km.labels_)

    km.fit(X)
    assert_array_equal(km.predict(X_csr), km.labels_)


@pytest.mark.parametrize("array_constr", [np.array, sp.csr_matrix],
                         ids=["dense", "sparse"])
@pytest.mark.parametrize("dtype", [np.int32, np.int64])
@pytest.mark.parametrize("init", ["k-means++", "ndarray"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_integer_input(Estimator, array_constr, dtype, init):
    # Check that KMeans and MiniBatchKMeans work with integer input.
    X_dense = np.array([[0, 0], [10, 10], [12, 9], [-1, 1], [2, 0], [8, 10]])
    X = array_constr(X_dense, dtype=dtype)

    n_init = 1 if init == "ndarray" else 10
    init = X_dense[:2] if init == "ndarray" else init

    km = Estimator(n_clusters=2, init=init, n_init=n_init, random_state=0)
    if Estimator is MiniBatchKMeans:
        km.set_params(batch_size=2)

    km.fit(X)

    # Internally integer input should be converted to float64
    assert km.cluster_centers_.dtype == np.float64

    expected_labels = [0, 1, 1, 0, 0, 1]
    assert_allclose(v_measure_score(km.labels_, expected_labels), 1)

    # Same with partial_fit (#14314)
    if Estimator is MiniBatchKMeans:
        km = clone(km).partial_fit(X)
        assert km.cluster_centers_.dtype == np.float64


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_transform(Estimator):
    # Check the transform method
    km = Estimator(n_clusters=n_clusters).fit(X)

    # Transorfming cluster_centers_ should return the pairwise distances
    # between centers
    Xt = km.transform(km.cluster_centers_)
    assert_allclose(Xt, pairwise_distances(km.cluster_centers_))
    # In particular, diagonal must be 0
    assert_array_equal(Xt.diagonal(), np.zeros(n_clusters))

    # Transorfming X should return the pairwise distances between X and the
    # centers
    Xt = km.transform(X)
    assert_allclose(Xt, pairwise_distances(X, km.cluster_centers_))


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_fit_transform(Estimator):
    # Check equivalence between fit.transform and fit_transform
    X1 = Estimator(random_state=0, n_init=1).fit(X).transform(X)
    X2 = Estimator(random_state=0, n_init=1).fit_transform(X)
    assert_allclose(X1, X2)


@pytest.mark.parametrize('algo', ['full', 'elkan'])
def test_predict_equal_labels(algo):
    km = KMeans(random_state=13, n_init=1, max_iter=1,
                algorithm=algo)
    km.fit(X)
    assert_array_equal(km.predict(X), km.labels_)


def test_n_init():
    # Check that increasing the number of init increases the quality
    previous_inertia = np.inf
    for n_init in [1, 5, 10]:
        # set max_iter=1 to avoid finding the global minimum and get the same
        # inertia each time
        km = KMeans(n_clusters=n_clusters, init="random", n_init=n_init,
                    random_state=0, max_iter=1).fit(X)
        assert km.inertia_ <= previous_inertia


def test_k_means_function():
    # test calling the k_means function directly
    cluster_centers, labels, inertia = k_means(X, n_clusters=n_clusters,
                                               sample_weight=None)

    assert cluster_centers.shape == (n_clusters, n_features)
    assert np.unique(labels).shape[0] == n_clusters

    # check that the labels assignment are perfect (up to a permutation)
    assert_allclose(v_measure_score(true_labels, labels), 1.0)
    assert inertia > 0.0


@pytest.mark.parametrize("data", [X, X_csr], ids=["dense", "sparse"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_float_precision(Estimator, data):
    # Check that the results are the same for single and double precision.
    km = Estimator(n_init=1, random_state=0)

    inertia = {}
    Xt = {}
    centers = {}
    labels = {}

    for dtype in [np.float64, np.float32]:
        X = data.astype(dtype, **_astype_copy_false(data))
        km.fit(X)

        inertia[dtype] = km.inertia_
        Xt[dtype] = km.transform(X)
        centers[dtype] = km.cluster_centers_
        labels[dtype] = km.labels_

        # dtype of cluster centers has to be the dtype of the input data
        assert km.cluster_centers_.dtype == dtype

        # same with partial_fit
        if Estimator is MiniBatchKMeans:
            km.partial_fit(X[0:3])
            assert km.cluster_centers_.dtype == dtype

    # compare arrays with low precision since the difference between 32 and
    # 64 bit comes from an accumulation of rounding errors.
    assert_allclose(inertia[np.float32], inertia[np.float64], rtol=1e-5)
    assert_allclose(Xt[np.float32], Xt[np.float64], rtol=1e-5)
    assert_allclose(centers[np.float32], centers[np.float64], rtol=1e-5)
    assert_array_equal(labels[np.float32], labels[np.float64])


@pytest.mark.parametrize("dtype", [np.int32, np.int64, np.float32, np.float64])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_centers_not_mutated(Estimator, dtype):
    # Check that KMeans and MiniBatchKMeans won't mutate the user provided
    # init centers silently even if input data and init centers have the same
    # type.
    X_new_type = X.astype(dtype, copy=False)
    centers_new_type = centers.astype(dtype, copy=False)

    km = Estimator(init=centers_new_type, n_clusters=n_clusters, n_init=1)
    km.fit(X_new_type)

    assert not np.may_share_memory(km.cluster_centers_, centers_new_type)


@pytest.mark.parametrize("data", [X, X_csr], ids=["dense", "sparse"])
def test_kmeans_init_fitted_centers(data):
    # Check that starting fitting from a local optimum shouldn't change the
    # solution
    km1 = KMeans(n_clusters=n_clusters).fit(data)
    km2 = KMeans(n_clusters=n_clusters, init=km1.cluster_centers_,
                 n_init=1).fit(data)

    assert_allclose(km1.cluster_centers_, km2.cluster_centers_)


def test_kmeans_warns_less_centers_than_unique_points():
    # Check KMeans when the number of found clusters is smaller than expected
    X = np.asarray([[0, 0],
                    [0, 1],
                    [1, 0],
                    [1, 0]])  # last point is duplicated
    km = KMeans(n_clusters=4)

    # KMeans should warn that fewer labels than cluster centers have been used
    msg = (r"Number of distinct clusters \(3\) found smaller than "
           r"n_clusters \(4\). Possibly due to duplicate points in X.")
    with pytest.warns(ConvergenceWarning, match=msg):
        km.fit(X)
        # only three distinct points, so only three clusters
        # can have points assigned to them
        assert set(km.labels_) == set(range(3))


def _sort_centers(centers):
    return np.sort(centers, axis=0)


def test_weighted_vs_repeated():
    # Check that a sample weight of N should yield the same result as an N-fold
    # repetition of the sample. Valid only if init is precomputed, otherwise
    # rng produces different results. Not valid for MinibatchKMeans due to rng
    # to extract minibatches.
    sample_weight = np.random.RandomState(0).randint(1, 5, size=n_samples)
    X_repeat = np.repeat(X, sample_weight, axis=0)

    km = KMeans(init=centers, n_init=1, n_clusters=n_clusters, random_state=0)

    km_weighted = clone(km).fit(X, sample_weight=sample_weight)
    repeated_labels = np.repeat(km_weighted.labels_, sample_weight)
    km_repeated = clone(km).fit(X_repeat)

    assert_array_equal(km_repeated.labels_, repeated_labels)
    assert_allclose(km_weighted.inertia_, km_repeated.inertia_)
    assert_allclose(_sort_centers(km_weighted.cluster_centers_),
                    _sort_centers(km_repeated.cluster_centers_))


@pytest.mark.parametrize("data", [X, X_csr], ids=["dense", "sparse"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_unit_weights_vs_no_weights(Estimator, data):
    # Check that not passing sample weights should be equivalent to passing
    # sample weights all equal to one.
    sample_weight = np.ones(n_samples)

    km = Estimator(n_clusters=n_clusters, random_state=42, n_init=1)
    km_none = clone(km).fit(data, sample_weight=None)
    km_ones = clone(km).fit(data, sample_weight=sample_weight)

    assert_array_equal(km_none.labels_, km_ones.labels_)
    assert_allclose(km_none.cluster_centers_, km_ones.cluster_centers_)


def test_scaled_weights():
    # scaling all sample weights by a common factor
    # shouldn't change the result
    sample_weight = np.ones(n_samples)
    for estimator in [KMeans(n_clusters=n_clusters, random_state=42),
                      MiniBatchKMeans(n_clusters=n_clusters, random_state=42)]:
        est_1 = clone(estimator).fit(X)
        est_2 = clone(estimator).fit(X, sample_weight=0.5*sample_weight)
        assert_almost_equal(v_measure_score(est_1.labels_, est_2.labels_), 1.0)
        assert_almost_equal(_sort_centers(est_1.cluster_centers_),
                            _sort_centers(est_2.cluster_centers_))


def test_kmeans_elkan_iter_attribute():
    # Regression test on bad n_iter_ value. Previous bug n_iter_ was one off
    # it's right value (#11340).
    km = KMeans(algorithm="elkan", max_iter=1).fit(X)
    assert km.n_iter_ == 1


def test_k_means_empty_cluster_relocated():
    # check that empty clusters are correctly relocated when using sample
    # weights (#13486)
    X = np.array([[-1], [1]])
    sample_weight = [1.9, 0.1]
    init = np.array([[-1], [10]])

    km = KMeans(n_clusters=2, init=init, n_init=1)
    km.fit(X, sample_weight=sample_weight)

    assert len(set(km.labels_)) == 2
    assert_allclose(km.cluster_centers_, [[-1], [1]])


def test_result_of_kmeans_equal_in_diff_n_threads():
    # Check that KMeans gives the same results in parallel mode than in
    # sequential mode.
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(50, 10))

    with threadpool_limits(limits=1, user_api="openmp"):
        result_1 = KMeans(
            n_clusters=3, random_state=0).fit(X).labels_
    with threadpool_limits(limits=2, user_api="openmp"):
        result_2 = KMeans(
            n_clusters=3, random_state=0).fit(X).labels_
    assert_array_equal(result_1, result_2)


@pytest.mark.parametrize("precompute_distances", ["auto", False, True])
def test_precompute_distance_deprecated(precompute_distances):
    # FIXME: remove in 0.25
    depr_msg = ("'precompute_distances' was deprecated in version 0.23 and "
                "will be removed in 0.25.")
    X, _ = make_blobs(n_samples=10, n_features=2, centers=2, random_state=0)
    kmeans = KMeans(n_clusters=2, n_init=1, init='random', random_state=0,
                    precompute_distances=precompute_distances)

    with pytest.warns(FutureWarning, match=depr_msg):
        kmeans.fit(X)


@pytest.mark.parametrize("n_jobs", [None, 1])
def test_n_jobs_deprecated(n_jobs):
    # FIXME: remove in 0.25
    depr_msg = ("'n_jobs' was deprecated in version 0.23 and will be removed "
                "in 0.25.")
    X, _ = make_blobs(n_samples=10, n_features=2, centers=2, random_state=0)
    kmeans = KMeans(n_clusters=2, n_init=1, init='random', random_state=0,
                    n_jobs=n_jobs)

    with pytest.warns(FutureWarning, match=depr_msg):
        kmeans.fit(X)


@pytest.mark.parametrize("attr", ["counts_", "init_size_", "random_state_"])
def test_minibatch_kmeans_deprecated_attributes(attr):
    # check that we raise a deprecation warning when accessing `init_size_`
    # FIXME: remove in 0.26
    depr_msg = (f"The attribute '{attr}' is deprecated in 0.24 and will be "
                f"removed in 0.26.")
    km = MiniBatchKMeans(n_clusters=2, n_init=1, init='random', random_state=0)
    km.fit(X)

    with pytest.warns(FutureWarning, match=depr_msg):
        getattr(km, attr)


def test_warning_elkan_1_cluster():
    X, _ = make_blobs(n_samples=10, n_features=2, centers=1, random_state=0)
    kmeans = KMeans(n_clusters=1, n_init=1, init='random', random_state=0,
                    algorithm='elkan')

    with pytest.warns(RuntimeWarning,
                      match="algorithm='elkan' doesn't make sense for a single"
                            " cluster"):
        kmeans.fit(X)


@pytest.mark.parametrize("array_constr",
                         [np.array, sp.csr_matrix],
                         ids=['dense', 'sparse'])
@pytest.mark.parametrize("algo", ['full', 'elkan'])
def test_k_means_1_iteration(array_constr, algo):
    # check the results after a single iteration (E-step M-step E-step) by
    # comparing against a pure python implementation.
    X = np.random.RandomState(0).uniform(size=(100, 5))
    init_centers = X[:5]
    X = array_constr(X)

    def py_kmeans(X, init):
        new_centers = init.copy()
        labels = pairwise_distances_argmin(X, init)
        for label in range(init.shape[0]):
            new_centers[label] = X[labels == label].mean(axis=0)
        labels = pairwise_distances_argmin(X, new_centers)
        return labels, new_centers

    py_labels, py_centers = py_kmeans(X, init_centers)

    cy_kmeans = KMeans(n_clusters=5, n_init=1, init=init_centers,
                       algorithm=algo, max_iter=1).fit(X)
    cy_labels = cy_kmeans.labels_
    cy_centers = cy_kmeans.cluster_centers_

    assert_array_equal(py_labels, cy_labels)
    assert_allclose(py_centers, cy_centers)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("squared", [True, False])
def test_euclidean_distance(dtype, squared):
    rng = np.random.RandomState(0)
    a_sparse = sp.random(1, 100, density=0.5, format="csr", random_state=rng,
                         dtype=dtype)
    a_dense = a_sparse.toarray().reshape(-1)
    b = rng.randn(100).astype(dtype, copy=False)
    b_squared_norm = (b**2).sum()

    expected = ((a_dense - b)**2).sum()
    expected = expected if squared else np.sqrt(expected)

    distance_dense_dense = _euclidean_dense_dense_wrapper(a_dense, b, squared)
    distance_sparse_dense = _euclidean_sparse_dense_wrapper(
        a_sparse.data, a_sparse.indices, b, b_squared_norm, squared)

    assert_allclose(distance_dense_dense, distance_sparse_dense, rtol=1e-6)
    assert_allclose(distance_dense_dense, expected, rtol=1e-6)
    assert_allclose(distance_sparse_dense, expected, rtol=1e-6)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_inertia(dtype):
    rng = np.random.RandomState(0)
    X_sparse = sp.random(100, 10, density=0.5, format="csr", random_state=rng,
                         dtype=dtype)
    X_dense = X_sparse.toarray()
    sample_weight = rng.randn(100).astype(dtype, copy=False)
    centers = rng.randn(5, 10).astype(dtype, copy=False)
    labels = rng.randint(5, size=100, dtype=np.int32)

    distances = ((X_dense - centers[labels])**2).sum(axis=1)
    expected = np.sum(distances * sample_weight)

    inertia_dense = _inertia_dense(X_dense, sample_weight, centers, labels)
    inertia_sparse = _inertia_sparse(X_sparse, sample_weight, centers, labels)

    assert_allclose(inertia_dense, inertia_sparse, rtol=1e-6)
    assert_allclose(inertia_dense, expected, rtol=1e-6)
    assert_allclose(inertia_sparse, expected, rtol=1e-6)


def test_sample_weight_unchanged():
    # Check that sample_weight is not modified in place by KMeans (#17204)
    X = np.array([[1], [2], [4]])
    sample_weight = np.array([0.5, 0.2, 0.3])
    KMeans(n_clusters=2, random_state=0).fit(X, sample_weight=sample_weight)

    assert_array_equal(sample_weight, np.array([0.5, 0.2, 0.3]))


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
@pytest.mark.parametrize("param, match", [
    ({"n_init": 0}, r"n_init should be > 0"),
    ({"max_iter": 0}, r"max_iter should be > 0"),
    ({"n_clusters": n_samples + 1}, r"n_samples.* should be >= n_clusters"),
    ({"init": X[:2]},
     r"The shape of the initial centers .* does not match "
     r"the number of clusters"),
    ({"init": lambda X_, k, random_state: X_[:2]},
     r"The shape of the initial centers .* does not match "
     r"the number of clusters"),
    ({"init": X[:8, :2]},
     r"The shape of the initial centers .* does not match "
     r"the number of features of the data"),
    ({"init": lambda X_, k, random_state: X_[:8, :2]},
     r"The shape of the initial centers .* does not match "
     r"the number of features of the data"),
    ({"init": "wrong"},
     r"init should be either 'k-means\+\+', 'random', "
     r"a ndarray or a callable")]
)
def test_wrong_params(Estimator, param, match):
    # Check that error are raised with clear error message when wrong values
    # are passed for the parameters
    with pytest.raises(ValueError, match=match):
        Estimator(**param).fit(X)


@pytest.mark.parametrize("param, match", [
    ({"algorithm": "wrong"}, r"Algorithm must be 'auto', 'full' or 'elkan'")]
)
def test_kmeans_wrong_params(param, match):
    # Check that error are raised with clear error message when wrong values
    # are passed for the KMeans specific parameters
    with pytest.raises(ValueError, match=match):
        KMeans(**param).fit(X)


@pytest.mark.parametrize("param, match", [
    ({"max_no_improvement": -1}, r"max_no_improvement should be >= 0"),
    ({"batch_size": -1}, r"batch_size should be > 0"),
    ({"init_size": -1}, r"init_size should be > 0"),
    ({"reassignment_ratio": -1}, r"reassignment_ratio should be >= 0")]
)
def test_minibatch_kmeans_wrong_params(param, match):
    # Check that error are raised with clear error message when wrong values
    # are passed for the MiniBatchKMeans specific parameters
    with pytest.raises(ValueError, match=match):
        MiniBatchKMeans(**param).fit(X)
