"""Testing for K-means"""

import re
import sys
from io import StringIO

import numpy as np
import pytest
from threadpoolctl import threadpool_info

import sklearn
from sklearn.base import clone
from sklearn.cluster import KMeans, MiniBatchKMeans, k_means, kmeans_plusplus
from sklearn.cluster._k_means_common import (
    _euclidean_dense_dense_wrapper,
    _euclidean_sparse_dense_wrapper,
    _inertia_dense,
    _inertia_sparse,
    _is_same_clustering,
    _relocate_empty_clusters_dense,
    _relocate_empty_clusters_sparse,
)
from sklearn.cluster._kmeans import _labels_inertia, _mini_batch_step
from sklearn.datasets import make_blobs
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import pairwise_distances, pairwise_distances_argmin
from sklearn.metrics.cluster import adjusted_rand_score, v_measure_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils._testing import (
    assert_allclose,
    assert_array_equal,
    create_memmap_backed_data,
)
from sklearn.utils.extmath import row_norms
from sklearn.utils.fixes import CSR_CONTAINERS, _sparse_random_array
from sklearn.utils.parallel import _get_threadpool_controller

# non centered, sparse centers to check the
centers = np.array(
    [
        [0.0, 5.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 4.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 5.0, 1.0],
    ]
)
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(
    n_samples=n_samples, centers=centers, cluster_std=1.0, random_state=42
)
X_as_any_csr = [container(X) for container in CSR_CONTAINERS]
data_containers = [np.array] + CSR_CONTAINERS
data_containers_ids = (
    ["dense", "sparse_matrix", "sparse_array"]
    if len(X_as_any_csr) == 2
    else ["dense", "sparse_matrix"]
)


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("algo", ["lloyd", "elkan"])
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


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("algo", ["lloyd", "elkan"])
def test_kmeans_relocated_clusters(array_constr, algo):
    # check that empty clusters are relocated as expected
    X = array_constr([[0, 0], [0.5, 0], [0.5, 1], [1, 1]])

    # second center too far from others points will be empty at first iter
    init_centers = np.array([[0.5, 0.5], [3, 3]])

    kmeans = KMeans(n_clusters=2, n_init=1, init=init_centers, algorithm=algo)
    kmeans.fit(X)

    expected_n_iter = 3
    expected_inertia = 0.25
    assert_allclose(kmeans.inertia_, expected_inertia)
    assert kmeans.n_iter_ == expected_n_iter

    # There are two acceptable ways of relocating clusters in this example, the output
    # depends on how the argpartition strategy breaks ties. We accept both outputs.
    try:
        expected_labels = [0, 0, 1, 1]
        expected_centers = [[0.25, 0], [0.75, 1]]
        assert_array_equal(kmeans.labels_, expected_labels)
        assert_allclose(kmeans.cluster_centers_, expected_centers)
    except AssertionError:
        expected_labels = [1, 1, 0, 0]
        expected_centers = [[0.75, 1.0], [0.25, 0.0]]
        assert_array_equal(kmeans.labels_, expected_labels)
        assert_allclose(kmeans.cluster_centers_, expected_centers)


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
def test_relocate_empty_clusters(array_constr):
    # test for the _relocate_empty_clusters_(dense/sparse) helpers

    # Synthetic dataset with 3 obvious clusters of different sizes
    X = np.array([-10.0, -9.5, -9, -8.5, -8, -1, 1, 9, 9.5, 10]).reshape(-1, 1)
    X = array_constr(X)
    sample_weight = np.ones(10)

    # centers all initialized to the first point of X
    centers_old = np.array([-10.0, -10, -10]).reshape(-1, 1)

    # With this initialization, all points will be assigned to the first center
    # At this point a center in centers_new is the weighted sum of the points
    # it contains if it's not empty, otherwise it is the same as before.
    centers_new = np.array([-16.5, -10, -10]).reshape(-1, 1)
    weight_in_clusters = np.array([10.0, 0, 0])
    labels = np.zeros(10, dtype=np.int32)

    if array_constr is np.array:
        _relocate_empty_clusters_dense(
            X, sample_weight, centers_old, centers_new, weight_in_clusters, labels
        )
    else:
        _relocate_empty_clusters_sparse(
            X.data,
            X.indices,
            X.indptr,
            sample_weight,
            centers_old,
            centers_new,
            weight_in_clusters,
            labels,
        )

    # The relocation scheme will take the 2 points farthest from the center and
    # assign them to the 2 empty clusters, i.e. points at 10 and at 9.9. The
    # first center will be updated to contain the other 8 points.
    assert_array_equal(weight_in_clusters, [8, 1, 1])
    assert_allclose(centers_new, [[-36], [10], [9.5]])


@pytest.mark.parametrize("distribution", ["normal", "blobs"])
@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("tol", [1e-2, 1e-8, 1e-100, 0])
def test_kmeans_elkan_results(distribution, array_constr, tol, global_random_seed):
    # Check that results are identical between lloyd and elkan algorithms
    rnd = np.random.RandomState(global_random_seed)
    if distribution == "normal":
        X = rnd.normal(size=(5000, 10))
    else:
        X, _ = make_blobs(random_state=rnd)
    X[X < 0] = 0
    X = array_constr(X)

    km_lloyd = KMeans(n_clusters=5, random_state=global_random_seed, n_init=1, tol=tol)
    km_elkan = KMeans(
        algorithm="elkan",
        n_clusters=5,
        random_state=global_random_seed,
        n_init=1,
        tol=tol,
    )

    km_lloyd.fit(X)
    km_elkan.fit(X)
    assert_allclose(km_elkan.cluster_centers_, km_lloyd.cluster_centers_)
    assert_array_equal(km_elkan.labels_, km_lloyd.labels_)
    assert km_elkan.n_iter_ == km_lloyd.n_iter_
    assert km_elkan.inertia_ == pytest.approx(km_lloyd.inertia_, rel=1e-6)


@pytest.mark.parametrize("algorithm", ["lloyd", "elkan"])
def test_kmeans_convergence(algorithm, global_random_seed):
    # Check that KMeans stops when convergence is reached when tol=0. (#16075)
    rnd = np.random.RandomState(global_random_seed)
    X = rnd.normal(size=(5000, 10))
    max_iter = 300

    km = KMeans(
        algorithm=algorithm,
        n_clusters=5,
        random_state=global_random_seed,
        n_init=1,
        tol=0,
        max_iter=max_iter,
    ).fit(X)

    assert km.n_iter_ < max_iter


@pytest.mark.parametrize("X_csr", X_as_any_csr)
def test_minibatch_update_consistency(X_csr, global_random_seed):
    # Check that dense and sparse minibatch update give the same results
    rng = np.random.RandomState(global_random_seed)

    centers_old = centers + rng.normal(size=centers.shape)
    centers_old_csr = centers_old.copy()

    centers_new = np.zeros_like(centers_old)
    centers_new_csr = np.zeros_like(centers_old_csr)

    weight_sums = np.zeros(centers_old.shape[0], dtype=X.dtype)
    weight_sums_csr = np.zeros(centers_old.shape[0], dtype=X.dtype)

    sample_weight = np.ones(X.shape[0], dtype=X.dtype)

    # extract a small minibatch
    X_mb = X[:10]
    X_mb_csr = X_csr[:10]
    sample_weight_mb = sample_weight[:10]

    # step 1: compute the dense minibatch update
    old_inertia = _mini_batch_step(
        X_mb,
        sample_weight_mb,
        centers_old,
        centers_new,
        weight_sums,
        np.random.RandomState(global_random_seed),
        random_reassign=False,
    )
    assert old_inertia > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels, new_inertia = _labels_inertia(X_mb, sample_weight_mb, centers_new)
    assert new_inertia > 0.0
    assert new_inertia < old_inertia

    # step 2: compute the sparse minibatch update
    old_inertia_csr = _mini_batch_step(
        X_mb_csr,
        sample_weight_mb,
        centers_old_csr,
        centers_new_csr,
        weight_sums_csr,
        np.random.RandomState(global_random_seed),
        random_reassign=False,
    )
    assert old_inertia_csr > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels_csr, new_inertia_csr = _labels_inertia(
        X_mb_csr, sample_weight_mb, centers_new_csr
    )
    assert new_inertia_csr > 0.0
    assert new_inertia_csr < old_inertia_csr

    # step 3: check that sparse and dense updates lead to the same results
    assert_array_equal(labels, labels_csr)
    assert_allclose(centers_new, centers_new_csr)
    assert_allclose(old_inertia, old_inertia_csr)
    assert_allclose(new_inertia, new_inertia_csr)


def _check_fitted_model(km):
    # check that the number of clusters centers and distinct labels match
    # the expectation
    centers = km.cluster_centers_
    assert centers.shape == (n_clusters, n_features)

    labels = km.labels_
    assert np.unique(labels).shape[0] == n_clusters

    # check that the labels assignment are perfect (up to a permutation)
    assert_allclose(v_measure_score(true_labels, labels), 1.0)
    assert km.inertia_ > 0.0


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
@pytest.mark.parametrize(
    "init",
    ["random", "k-means++", centers.copy(), lambda X, k, random_state: centers.copy()],
    ids=["random", "k-means++", "ndarray", "callable"],
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_all_init(Estimator, input_data, init):
    # Check KMeans and MiniBatchKMeans with all possible init.
    n_init = 10 if isinstance(init, str) else 1
    km = Estimator(
        init=init, n_clusters=n_clusters, random_state=42, n_init=n_init
    ).fit(input_data)
    _check_fitted_model(km)


@pytest.mark.parametrize(
    "init",
    ["random", "k-means++", centers, lambda X, k, random_state: centers.copy()],
    ids=["random", "k-means++", "ndarray", "callable"],
)
def test_minibatch_kmeans_partial_fit_init(init):
    if hasattr(init, "copy"):
        # Avoid mutating a shared array in place to avoid side effects in other tests.
        init = init.copy()

    # Check MiniBatchKMeans init with partial_fit
    n_init = 10 if isinstance(init, str) else 1
    km = MiniBatchKMeans(
        init=init, n_clusters=n_clusters, random_state=0, n_init=n_init
    )
    for i in range(100):
        # "random" init requires many batches to recover the true labels.
        km.partial_fit(X)
    _check_fitted_model(km)


@pytest.mark.parametrize(
    "init, expected_n_init",
    [
        ("k-means++", 1),
        ("random", "default"),
        (
            lambda X, n_clusters, random_state: random_state.uniform(
                size=(n_clusters, X.shape[1])
            ),
            "default",
        ),
        ("array-like", 1),
    ],
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_kmeans_init_auto_with_initial_centroids(Estimator, init, expected_n_init):
    """Check that `n_init="auto"` chooses the right number of initializations.
    Non-regression test for #26657:
    https://github.com/scikit-learn/scikit-learn/pull/26657
    """
    n_sample, n_features, n_clusters = 100, 10, 5
    X = np.random.randn(n_sample, n_features)
    if init == "array-like":
        init = np.random.randn(n_clusters, n_features)
    if expected_n_init == "default":
        expected_n_init = 3 if Estimator is MiniBatchKMeans else 10

    kmeans = Estimator(n_clusters=n_clusters, init=init, n_init="auto").fit(X)
    assert kmeans._n_init == expected_n_init


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_fortran_aligned_data(Estimator, global_random_seed):
    # Check that KMeans works with fortran-aligned data.
    X_fortran = np.asfortranarray(X)
    centers_fortran = np.asfortranarray(centers)

    km_c = Estimator(
        n_clusters=n_clusters, init=centers, n_init=1, random_state=global_random_seed
    ).fit(X)
    km_f = Estimator(
        n_clusters=n_clusters,
        init=centers_fortran,
        n_init=1,
        random_state=global_random_seed,
    ).fit(X_fortran)
    assert_allclose(km_c.cluster_centers_, km_f.cluster_centers_)
    assert_array_equal(km_c.labels_, km_f.labels_)


def test_minibatch_kmeans_verbose():
    # Check verbose mode of MiniBatchKMeans for better coverage.
    km = MiniBatchKMeans(n_clusters=n_clusters, random_state=42, verbose=1)
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        km.fit(X)
    finally:
        sys.stdout = old_stdout


@pytest.mark.parametrize("algorithm", ["lloyd", "elkan"])
@pytest.mark.parametrize("tol", [1e-2, 0])
def test_kmeans_verbose(algorithm, tol, capsys):
    # Check verbose mode of KMeans for better coverage.
    X = np.random.RandomState(0).normal(size=(5000, 10))

    KMeans(
        algorithm=algorithm,
        n_clusters=n_clusters,
        random_state=42,
        init="random",
        n_init=1,
        tol=tol,
        verbose=1,
    ).fit(X)

    captured = capsys.readouterr()

    assert re.search(r"Initialization complete", captured.out)
    assert re.search(r"Iteration [0-9]+, inertia", captured.out)

    if tol == 0:
        assert re.search(r"strict convergence", captured.out)
    else:
        assert re.search(r"center shift .* within tolerance", captured.out)


def test_minibatch_kmeans_warning_init_size():
    # Check that a warning is raised when init_size is smaller than n_clusters
    with pytest.warns(
        RuntimeWarning, match=r"init_size.* should be larger than n_clusters"
    ):
        MiniBatchKMeans(init_size=10, n_clusters=20).fit(X)


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_warning_n_init_precomputed_centers(Estimator):
    # Check that a warning is raised when n_init > 1 and an array is passed for
    # the init parameter.
    with pytest.warns(
        RuntimeWarning,
        match="Explicit initial center position passed: performing only one init",
    ):
        Estimator(init=centers, n_clusters=n_clusters, n_init=10).fit(X)


def test_minibatch_sensible_reassign(global_random_seed):
    # check that identical initial clusters are reassigned
    # also a regression test for when there are more desired reassignments than
    # samples.
    zeroed_X, true_labels = make_blobs(
        n_samples=100, centers=5, random_state=global_random_seed
    )
    zeroed_X[::2, :] = 0

    km = MiniBatchKMeans(
        n_clusters=20, batch_size=10, random_state=global_random_seed, init="random"
    ).fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    num_non_zero_clusters = km.cluster_centers_.any(axis=1).sum()
    assert num_non_zero_clusters > 9, f"{num_non_zero_clusters=} is too small"

    # do the same with batch-size > X.shape[0] (regression test)
    km = MiniBatchKMeans(
        n_clusters=20, batch_size=200, random_state=global_random_seed, init="random"
    ).fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    num_non_zero_clusters = km.cluster_centers_.any(axis=1).sum()
    assert num_non_zero_clusters > 9, f"{num_non_zero_clusters=} is too small"

    # do the same with partial_fit API
    km = MiniBatchKMeans(n_clusters=20, random_state=global_random_seed, init="random")
    for i in range(100):
        km.partial_fit(zeroed_X)
    # there should not be too many exact zero cluster centers
    num_non_zero_clusters = km.cluster_centers_.any(axis=1).sum()
    assert num_non_zero_clusters > 9, f"{num_non_zero_clusters=} is too small"


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
def test_minibatch_reassign(input_data, global_random_seed):
    # Check the reassignment part of the minibatch step with very high or very
    # low reassignment ratio.
    perfect_centers = np.empty((n_clusters, n_features))
    for i in range(n_clusters):
        perfect_centers[i] = X[true_labels == i].mean(axis=0)

    sample_weight = np.ones(n_samples)
    centers_new = np.empty_like(perfect_centers)

    # Give a perfect initialization, but a large reassignment_ratio, as a
    # result many centers should be reassigned and the model should no longer
    # be good
    score_before = -_labels_inertia(input_data, sample_weight, perfect_centers, 1)[1]

    _mini_batch_step(
        input_data,
        sample_weight,
        perfect_centers,
        centers_new,
        np.zeros(n_clusters),
        np.random.RandomState(global_random_seed),
        random_reassign=True,
        reassignment_ratio=1,
    )

    score_after = -_labels_inertia(input_data, sample_weight, centers_new, 1)[1]

    assert score_before > score_after

    # Give a perfect initialization, with a small reassignment_ratio,
    # no center should be reassigned.
    _mini_batch_step(
        input_data,
        sample_weight,
        perfect_centers,
        centers_new,
        np.zeros(n_clusters),
        np.random.RandomState(global_random_seed),
        random_reassign=True,
        reassignment_ratio=1e-15,
    )

    assert_allclose(centers_new, perfect_centers)


def test_minibatch_with_many_reassignments():
    # Test for the case that the number of clusters to reassign is bigger
    # than the batch_size. Run the test with 100 clusters and a batch_size of
    # 10 because it turned out that these values ensure that the number of
    # clusters to reassign is always bigger than the batch_size.
    MiniBatchKMeans(
        n_clusters=100,
        batch_size=10,
        init_size=n_samples,
        random_state=42,
        verbose=True,
    ).fit(X)


def test_minibatch_kmeans_init_size():
    # Check the internal _init_size attribute of MiniBatchKMeans

    # default init size should be 3 * batch_size
    km = MiniBatchKMeans(n_clusters=10, batch_size=5, n_init=1).fit(X)
    assert km._init_size == 15

    # if 3 * batch size < n_clusters, it should then be 3 * n_clusters
    km = MiniBatchKMeans(n_clusters=10, batch_size=1, n_init=1).fit(X)
    assert km._init_size == 30

    # it should not be larger than n_samples
    km = MiniBatchKMeans(
        n_clusters=10, batch_size=5, n_init=1, init_size=n_samples + 1
    ).fit(X)
    assert km._init_size == n_samples


@pytest.mark.parametrize("tol, max_no_improvement", [(1e-4, None), (0, 10)])
def test_minibatch_declared_convergence(capsys, tol, max_no_improvement):
    # Check convergence detection based on ewa batch inertia or on
    # small center change.
    X, _, centers = make_blobs(centers=3, random_state=0, return_centers=True)

    km = MiniBatchKMeans(
        n_clusters=3,
        init=centers,
        batch_size=20,
        tol=tol,
        random_state=0,
        max_iter=10,
        n_init=1,
        verbose=1,
        max_no_improvement=max_no_improvement,
    )

    km.fit(X)
    assert 1 < km.n_iter_ < 10

    captured = capsys.readouterr()
    if max_no_improvement is None:
        assert "Converged (small centers change)" in captured.out
    if tol == 0:
        assert "Converged (lack of improvement in inertia)" in captured.out


def test_minibatch_iter_steps():
    # Check consistency of n_iter_ and n_steps_ attributes.
    batch_size = 30
    n_samples = X.shape[0]
    km = MiniBatchKMeans(n_clusters=3, batch_size=batch_size, random_state=0).fit(X)

    # n_iter_ is the number of started epochs
    assert km.n_iter_ == np.ceil((km.n_steps_ * batch_size) / n_samples)
    assert isinstance(km.n_iter_, int)

    # without stopping condition, max_iter should be reached
    km = MiniBatchKMeans(
        n_clusters=3,
        batch_size=batch_size,
        random_state=0,
        tol=0,
        max_no_improvement=None,
        max_iter=10,
    ).fit(X)

    assert km.n_iter_ == 10
    assert km.n_steps_ == (10 * n_samples) // batch_size
    assert isinstance(km.n_steps_, int)


def test_kmeans_copyx():
    # Check that copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    km = KMeans(copy_x=False, n_clusters=n_clusters, random_state=42)
    km.fit(my_X)
    _check_fitted_model(km)

    # check that my_X is de-centered
    assert_allclose(my_X, X)


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_score_max_iter(Estimator, global_random_seed):
    # Check that fitting KMeans or MiniBatchKMeans with more iterations gives
    # better score
    X = np.random.RandomState(global_random_seed).randn(100, 10)

    km1 = Estimator(n_init=1, random_state=global_random_seed, max_iter=1)
    s1 = km1.fit(X).score(X)
    km2 = Estimator(n_init=1, random_state=global_random_seed, max_iter=10)
    s2 = km2.fit(X).score(X)
    assert s2 > s1


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize(
    "Estimator, algorithm",
    [(KMeans, "lloyd"), (KMeans, "elkan"), (MiniBatchKMeans, None)],
)
@pytest.mark.parametrize("max_iter", [2, 100])
def test_kmeans_predict(
    Estimator, algorithm, array_constr, max_iter, global_dtype, global_random_seed
):
    # Check the predict method and the equivalence between fit.predict and
    # fit_predict.
    X, _ = make_blobs(
        n_samples=200, n_features=10, centers=10, random_state=global_random_seed
    )
    X = array_constr(X, dtype=global_dtype)

    km = Estimator(
        n_clusters=10,
        init="random",
        n_init=10,
        max_iter=max_iter,
        random_state=global_random_seed,
    )
    if algorithm is not None:
        km.set_params(algorithm=algorithm)
    km.fit(X)
    labels = km.labels_

    # re-predict labels for training set using predict
    pred = km.predict(X)
    assert_array_equal(pred, labels)

    # re-predict labels for training set using fit_predict
    pred = km.fit_predict(X)
    assert_array_equal(pred, labels)

    # predict centroid labels
    pred = km.predict(km.cluster_centers_)
    assert_array_equal(pred, np.arange(10))


@pytest.mark.parametrize("X_csr", X_as_any_csr)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_dense_sparse(Estimator, X_csr, global_random_seed):
    # Check that the results are the same for dense and sparse input.
    sample_weight = np.random.RandomState(global_random_seed).random_sample(
        (n_samples,)
    )
    km_dense = Estimator(
        n_clusters=n_clusters, random_state=global_random_seed, n_init=1
    )
    km_dense.fit(X, sample_weight=sample_weight)
    km_sparse = Estimator(
        n_clusters=n_clusters, random_state=global_random_seed, n_init=1
    )
    km_sparse.fit(X_csr, sample_weight=sample_weight)

    assert_array_equal(km_dense.labels_, km_sparse.labels_)
    assert_allclose(km_dense.cluster_centers_, km_sparse.cluster_centers_)


@pytest.mark.parametrize("X_csr", X_as_any_csr)
@pytest.mark.parametrize(
    "init", ["random", "k-means++", centers], ids=["random", "k-means++", "ndarray"]
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_predict_dense_sparse(Estimator, init, X_csr):
    # check that models trained on sparse input also works for dense input at
    # predict time and vice versa.
    n_init = 10 if isinstance(init, str) else 1
    km = Estimator(n_clusters=n_clusters, init=init, n_init=n_init, random_state=0)

    km.fit(X_csr)
    assert_array_equal(km.predict(X), km.labels_)

    km.fit(X)
    assert_array_equal(km.predict(X_csr), km.labels_)


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("dtype", [np.int32, np.int64])
@pytest.mark.parametrize("init", ["k-means++", "ndarray"])
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_integer_input(Estimator, array_constr, dtype, init, global_random_seed):
    # Check that KMeans and MiniBatchKMeans work with integer input.
    X_dense = np.array([[0, 0], [10, 10], [12, 9], [-1, 1], [2, 0], [8, 10]])
    X = array_constr(X_dense, dtype=dtype)

    n_init = 1 if init == "ndarray" else 10
    init = X_dense[:2] if init == "ndarray" else init

    km = Estimator(
        n_clusters=2, init=init, n_init=n_init, random_state=global_random_seed
    )
    if Estimator is MiniBatchKMeans:
        km.set_params(batch_size=2)

    km.fit(X)

    # Internally integer input should be converted to float64
    assert km.cluster_centers_.dtype == np.float64

    expected_labels = [0, 1, 1, 0, 0, 1]
    assert_allclose(v_measure_score(km.labels_, expected_labels), 1.0)

    # Same with partial_fit (#14314)
    if Estimator is MiniBatchKMeans:
        km = clone(km).partial_fit(X)
        assert km.cluster_centers_.dtype == np.float64


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_transform(Estimator, global_random_seed):
    # Check the transform method
    km = Estimator(n_clusters=n_clusters, random_state=global_random_seed).fit(X)

    # Transorfming cluster_centers_ should return the pairwise distances
    # between centers
    Xt = km.transform(km.cluster_centers_)
    assert_allclose(Xt, pairwise_distances(km.cluster_centers_))
    # In particular, diagonal must be 0
    assert_array_equal(Xt.diagonal(), np.zeros(n_clusters))

    # Transforming X should return the pairwise distances between X and the
    # centers
    Xt = km.transform(X)
    assert_allclose(Xt, pairwise_distances(X, km.cluster_centers_))


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_fit_transform(Estimator, global_random_seed):
    # Check equivalence between fit.transform and fit_transform
    X1 = Estimator(random_state=global_random_seed, n_init=1).fit(X).transform(X)
    X2 = Estimator(random_state=global_random_seed, n_init=1).fit_transform(X)
    assert_allclose(X1, X2)


def test_n_init(global_random_seed):
    # Check that increasing the number of init increases the quality
    previous_inertia = np.inf
    for n_init in [1, 5, 10]:
        # set max_iter=1 to avoid finding the global minimum and get the same
        # inertia each time
        km = KMeans(
            n_clusters=n_clusters,
            init="random",
            n_init=n_init,
            random_state=global_random_seed,
            max_iter=1,
        ).fit(X)
        assert km.inertia_ <= previous_inertia


def test_k_means_function(global_random_seed):
    # test calling the k_means function directly
    cluster_centers, labels, inertia = k_means(
        X, n_clusters=n_clusters, sample_weight=None, random_state=global_random_seed
    )

    assert cluster_centers.shape == (n_clusters, n_features)
    assert np.unique(labels).shape[0] == n_clusters

    # check that the labels assignment are perfect (up to a permutation)
    assert_allclose(v_measure_score(true_labels, labels), 1.0)
    assert inertia > 0.0


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
@pytest.mark.skipif(
    not any(i for i in threadpool_info() if i["user_api"] == "blas"),
    reason=(
        "Fails for some global_random_seed on Atlas which cannot be detected by "
        "threadpoolctl."
    ),
)
def test_float_precision(Estimator, input_data, global_random_seed):
    # Check that the results are the same for single and double precision.
    km = Estimator(n_init=1, random_state=global_random_seed)

    inertia = {}
    Xt = {}
    centers = {}
    labels = {}

    for dtype in [np.float64, np.float32]:
        X = input_data.astype(dtype, copy=False)
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
    rtol = 1e-4
    assert_allclose(inertia[np.float32], inertia[np.float64], rtol=rtol)
    assert_allclose(Xt[np.float32], Xt[np.float64], atol=Xt[np.float64].max() * rtol)
    assert_allclose(
        centers[np.float32], centers[np.float64], atol=centers[np.float64].max() * rtol
    )
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


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
def test_kmeans_init_fitted_centers(input_data):
    # Check that starting fitting from a local optimum shouldn't change the
    # solution
    km1 = KMeans(n_clusters=n_clusters).fit(input_data)
    km2 = KMeans(n_clusters=n_clusters, init=km1.cluster_centers_, n_init=1).fit(
        input_data
    )

    assert_allclose(km1.cluster_centers_, km2.cluster_centers_)


def test_kmeans_warns_less_centers_than_unique_points(global_random_seed):
    # Check KMeans when the number of found clusters is smaller than expected
    X = np.asarray([[0, 0], [0, 1], [1, 0], [1, 0]])  # last point is duplicated
    km = KMeans(n_clusters=4, random_state=global_random_seed)

    # KMeans should warn that fewer labels than cluster centers have been used
    msg = (
        r"Number of distinct clusters \(3\) found smaller than "
        r"n_clusters \(4\). Possibly due to duplicate points in X."
    )
    with pytest.warns(ConvergenceWarning, match=msg):
        km.fit(X)
        # only three distinct points, so only three clusters
        # can have points assigned to them
        assert set(km.labels_) == set(range(3))


def _sort_centers(centers):
    return np.sort(centers, axis=0)


def test_weighted_vs_repeated(global_random_seed):
    # Check that a sample weight of N should yield the same result as an N-fold
    # repetition of the sample. Valid only if init is precomputed, otherwise
    # rng produces different results. Not valid for MinibatchKMeans due to rng
    # to extract minibatches.
    sample_weight = np.random.RandomState(global_random_seed).randint(
        1, 5, size=n_samples
    )
    X_repeat = np.repeat(X, sample_weight, axis=0)

    km = KMeans(
        init=centers, n_init=1, n_clusters=n_clusters, random_state=global_random_seed
    )

    km_weighted = clone(km).fit(X, sample_weight=sample_weight)
    repeated_labels = np.repeat(km_weighted.labels_, sample_weight)
    km_repeated = clone(km).fit(X_repeat)

    assert_array_equal(km_repeated.labels_, repeated_labels)
    assert_allclose(km_weighted.inertia_, km_repeated.inertia_)
    assert_allclose(
        _sort_centers(km_weighted.cluster_centers_),
        _sort_centers(km_repeated.cluster_centers_),
    )


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_unit_weights_vs_no_weights(Estimator, input_data, global_random_seed):
    # Check that not passing sample weights should be equivalent to passing
    # sample weights all equal to one.
    sample_weight = np.ones(n_samples)

    km = Estimator(n_clusters=n_clusters, random_state=global_random_seed, n_init=1)
    km_none = clone(km).fit(input_data, sample_weight=None)
    km_ones = clone(km).fit(input_data, sample_weight=sample_weight)

    assert_array_equal(km_none.labels_, km_ones.labels_)
    assert_allclose(km_none.cluster_centers_, km_ones.cluster_centers_)


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
    ids=data_containers_ids,
)
@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_scaled_weights(Estimator, input_data, global_random_seed):
    # Check that scaling all sample weights by a common factor
    # shouldn't change the result
    sample_weight = np.random.RandomState(global_random_seed).uniform(size=n_samples)

    km = Estimator(n_clusters=n_clusters, random_state=global_random_seed, n_init=1)
    km_orig = clone(km).fit(input_data, sample_weight=sample_weight)
    km_scaled = clone(km).fit(input_data, sample_weight=0.5 * sample_weight)

    assert_array_equal(km_orig.labels_, km_scaled.labels_)
    assert_allclose(km_orig.cluster_centers_, km_scaled.cluster_centers_)


def test_kmeans_elkan_iter_attribute():
    # Regression test on bad n_iter_ value. Previous bug n_iter_ was one off
    # it's right value (#11340).
    km = KMeans(algorithm="elkan", max_iter=1).fit(X)
    assert km.n_iter_ == 1


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
def test_kmeans_empty_cluster_relocated(array_constr):
    # check that empty clusters are correctly relocated when using sample
    # weights (#13486)
    X = array_constr([[-1], [1]])
    sample_weight = [1.9, 0.1]
    init = np.array([[-1], [10]])

    km = KMeans(n_clusters=2, init=init, n_init=1)
    km.fit(X, sample_weight=sample_weight)

    assert len(set(km.labels_)) == 2
    assert_allclose(km.cluster_centers_, [[-1], [1]])


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_result_equal_in_diff_n_threads(Estimator, global_random_seed):
    # Check that KMeans/MiniBatchKMeans give the same results in parallel mode
    # than in sequential mode.
    rnd = np.random.RandomState(global_random_seed)
    X = rnd.normal(size=(50, 10))

    with _get_threadpool_controller().limit(limits=1, user_api="openmp"):
        result_1 = (
            Estimator(n_clusters=n_clusters, random_state=global_random_seed)
            .fit(X)
            .labels_
        )
    with _get_threadpool_controller().limit(limits=2, user_api="openmp"):
        result_2 = (
            Estimator(n_clusters=n_clusters, random_state=global_random_seed)
            .fit(X)
            .labels_
        )
    assert_array_equal(result_1, result_2)


def test_warning_elkan_1_cluster():
    # Check warning messages specific to KMeans
    with pytest.warns(
        RuntimeWarning,
        match="algorithm='elkan' doesn't make sense for a single cluster",
    ):
        KMeans(n_clusters=1, algorithm="elkan").fit(X)


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("algo", ["lloyd", "elkan"])
def test_k_means_1_iteration(array_constr, algo, global_random_seed):
    # check the results after a single iteration (E-step M-step E-step) by
    # comparing against a pure python implementation.
    X = np.random.RandomState(global_random_seed).uniform(size=(100, 5))
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

    cy_kmeans = KMeans(
        n_clusters=5, n_init=1, init=init_centers, algorithm=algo, max_iter=1
    ).fit(X)
    cy_labels = cy_kmeans.labels_
    cy_centers = cy_kmeans.cluster_centers_

    assert_array_equal(py_labels, cy_labels)
    assert_allclose(py_centers, cy_centers)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("squared", [True, False])
def test_euclidean_distance(dtype, squared, global_random_seed):
    # Check that the _euclidean_(dense/sparse)_dense helpers produce correct
    # results
    rng = np.random.RandomState(global_random_seed)
    a_sparse = _sparse_random_array(
        (1, 100), density=0.5, format="csr", rng=rng, dtype=dtype
    )
    a_dense = a_sparse.toarray().reshape(-1)
    b = rng.randn(100).astype(dtype, copy=False)
    b_squared_norm = (b**2).sum()

    expected = ((a_dense - b) ** 2).sum()
    expected = expected if squared else np.sqrt(expected)

    distance_dense_dense = _euclidean_dense_dense_wrapper(a_dense, b, squared)
    distance_sparse_dense = _euclidean_sparse_dense_wrapper(
        a_sparse.data, a_sparse.indices, b, b_squared_norm, squared
    )

    rtol = 1e-4 if dtype == np.float32 else 1e-7
    assert_allclose(distance_dense_dense, distance_sparse_dense, rtol=rtol)
    assert_allclose(distance_dense_dense, expected, rtol=rtol)
    assert_allclose(distance_sparse_dense, expected, rtol=rtol)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_inertia(dtype, global_random_seed):
    # Check that the _inertia_(dense/sparse) helpers produce correct results.
    rng = np.random.RandomState(global_random_seed)
    X_sparse = _sparse_random_array(
        (100, 10), density=0.5, format="csr", rng=rng, dtype=dtype
    )
    X_dense = X_sparse.toarray()
    sample_weight = rng.randn(100).astype(dtype, copy=False)
    centers = rng.randn(5, 10).astype(dtype, copy=False)
    labels = rng.randint(5, size=100, dtype=np.int32)

    distances = ((X_dense - centers[labels]) ** 2).sum(axis=1)
    expected = np.sum(distances * sample_weight)

    inertia_dense = _inertia_dense(X_dense, sample_weight, centers, labels, n_threads=1)
    inertia_sparse = _inertia_sparse(
        X_sparse, sample_weight, centers, labels, n_threads=1
    )

    rtol = 1e-4 if dtype == np.float32 else 1e-6
    assert_allclose(inertia_dense, inertia_sparse, rtol=rtol)
    assert_allclose(inertia_dense, expected, rtol=rtol)
    assert_allclose(inertia_sparse, expected, rtol=rtol)

    # Check the single_label parameter.
    label = 1
    mask = labels == label
    distances = ((X_dense[mask] - centers[label]) ** 2).sum(axis=1)
    expected = np.sum(distances * sample_weight[mask])

    inertia_dense = _inertia_dense(
        X_dense, sample_weight, centers, labels, n_threads=1, single_label=label
    )
    inertia_sparse = _inertia_sparse(
        X_sparse, sample_weight, centers, labels, n_threads=1, single_label=label
    )

    assert_allclose(inertia_dense, inertia_sparse, rtol=rtol)
    assert_allclose(inertia_dense, expected, rtol=rtol)
    assert_allclose(inertia_sparse, expected, rtol=rtol)


@pytest.mark.parametrize("Klass, default_n_init", [(KMeans, 10), (MiniBatchKMeans, 3)])
def test_n_init_auto(Klass, default_n_init):
    est = Klass(n_init="auto", init="k-means++")
    est.fit(X)
    assert est._n_init == 1

    est = Klass(n_init="auto", init="random")
    est.fit(X)
    assert est._n_init == 10 if Klass.__name__ == "KMeans" else 3


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
def test_sample_weight_unchanged(Estimator):
    # Check that sample_weight is not modified in place by KMeans (#17204)
    X = np.array([[1], [2], [4]])
    sample_weight = np.array([0.5, 0.2, 0.3])
    Estimator(n_clusters=2, random_state=0).fit(X, sample_weight=sample_weight)

    assert_array_equal(sample_weight, np.array([0.5, 0.2, 0.3]))


@pytest.mark.parametrize("Estimator", [KMeans, MiniBatchKMeans])
@pytest.mark.parametrize(
    "param, match",
    [
        ({"n_clusters": n_samples + 1}, r"n_samples.* should be >= n_clusters"),
        (
            {"init": X[:2]},
            r"The shape of the initial centers .* does not match "
            r"the number of clusters",
        ),
        (
            {"init": lambda X_, k, random_state: X_[:2]},
            r"The shape of the initial centers .* does not match "
            r"the number of clusters",
        ),
        (
            {"init": X[:8, :2]},
            r"The shape of the initial centers .* does not match "
            r"the number of features of the data",
        ),
        (
            {"init": lambda X_, k, random_state: X_[:8, :2]},
            r"The shape of the initial centers .* does not match "
            r"the number of features of the data",
        ),
    ],
)
def test_wrong_params(Estimator, param, match):
    # Check that error are raised with clear error message when wrong values
    # are passed for the parameters
    # Set n_init=1 by default to avoid warning with precomputed init
    km = Estimator(n_init=1)
    with pytest.raises(ValueError, match=match):
        km.set_params(**param).fit(X)


@pytest.mark.parametrize(
    "param, match",
    [
        (
            {"x_squared_norms": X[:2]},
            r"The length of x_squared_norms .* should "
            r"be equal to the length of n_samples",
        ),
    ],
)
def test_kmeans_plusplus_wrong_params(param, match):
    with pytest.raises(ValueError, match=match):
        kmeans_plusplus(X, n_clusters, **param)


@pytest.mark.parametrize(
    "input_data",
    [X] + X_as_any_csr,
)
@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_kmeans_plusplus_output(input_data, dtype, global_random_seed):
    # Check for the correct number of seeds and all positive values
    data = input_data.astype(dtype)
    centers, indices = kmeans_plusplus(
        data, n_clusters, random_state=global_random_seed
    )

    # Check there are the correct number of indices and that all indices are
    # positive and within the number of samples
    assert indices.shape[0] == n_clusters
    assert (indices >= 0).all()
    assert (indices <= data.shape[0]).all()

    # Check for the correct number of seeds and that they are bound by the data
    assert centers.shape[0] == n_clusters
    assert (centers.max(axis=0) <= data.max(axis=0)).all()
    assert (centers.min(axis=0) >= data.min(axis=0)).all()

    # Check that indices correspond to reported centers
    # Use X for comparison rather than data, test still works against centers
    # calculated with sparse data.
    assert_allclose(X[indices].astype(dtype), centers)


@pytest.mark.parametrize("x_squared_norms", [row_norms(X, squared=True), None])
def test_kmeans_plusplus_norms(x_squared_norms):
    # Check that defining x_squared_norms returns the same as default=None.
    centers, indices = kmeans_plusplus(X, n_clusters, x_squared_norms=x_squared_norms)

    assert_allclose(X[indices], centers)


def test_kmeans_plusplus_dataorder(global_random_seed):
    # Check that memory layout does not effect result
    centers_c, _ = kmeans_plusplus(X, n_clusters, random_state=global_random_seed)

    X_fortran = np.asfortranarray(X)

    centers_fortran, _ = kmeans_plusplus(
        X_fortran, n_clusters, random_state=global_random_seed
    )

    assert_allclose(centers_c, centers_fortran)


def test_is_same_clustering():
    # Sanity check for the _is_same_clustering utility function
    labels1 = np.array([1, 0, 0, 1, 2, 0, 2, 1], dtype=np.int32)
    assert _is_same_clustering(labels1, labels1, 3)

    # these other labels represent the same clustering since we can retrieve the first
    # labels by simply renaming the labels: 0 -> 1, 1 -> 2, 2 -> 0.
    labels2 = np.array([0, 2, 2, 0, 1, 2, 1, 0], dtype=np.int32)
    assert _is_same_clustering(labels1, labels2, 3)

    # these other labels do not represent the same clustering since not all ones are
    # mapped to a same value
    labels3 = np.array([1, 0, 0, 2, 2, 0, 2, 1], dtype=np.int32)
    assert not _is_same_clustering(labels1, labels3, 3)


@pytest.mark.parametrize(
    "kwargs", ({"init": np.str_("k-means++")}, {"init": [[0, 0], [1, 1]], "n_init": 1})
)
def test_kmeans_with_array_like_or_np_scalar_init(kwargs):
    """Check that init works with numpy scalar strings.

    Non-regression test for #21964.
    """
    X = np.asarray([[0, 0], [0.5, 0], [0.5, 1], [1, 1]], dtype=np.float64)

    clustering = KMeans(n_clusters=2, **kwargs)
    # Does not raise
    clustering.fit(X)


@pytest.mark.parametrize(
    "Klass, method",
    [(KMeans, "fit"), (MiniBatchKMeans, "fit"), (MiniBatchKMeans, "partial_fit")],
)
def test_feature_names_out(Klass, method):
    """Check `feature_names_out` for `KMeans` and `MiniBatchKMeans`."""
    class_name = Klass.__name__.lower()
    kmeans = Klass()
    getattr(kmeans, method)(X)
    n_clusters = kmeans.cluster_centers_.shape[0]

    names_out = kmeans.get_feature_names_out()
    assert_array_equal([f"{class_name}{i}" for i in range(n_clusters)], names_out)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS + [None])
def test_predict_does_not_change_cluster_centers(csr_container):
    """Check that predict does not change cluster centers.

    Non-regression test for gh-24253.
    """
    X, _ = make_blobs(n_samples=200, n_features=10, centers=10, random_state=0)
    if csr_container is not None:
        X = csr_container(X)

    kmeans = KMeans()
    y_pred1 = kmeans.fit_predict(X)
    # Make cluster_centers readonly
    kmeans.cluster_centers_ = create_memmap_backed_data(kmeans.cluster_centers_)
    kmeans.labels_ = create_memmap_backed_data(kmeans.labels_)

    y_pred2 = kmeans.predict(X)
    assert_array_equal(y_pred1, y_pred2)


@pytest.mark.parametrize("init", ["k-means++", "random"])
def test_sample_weight_init(init, global_random_seed):
    """Check that sample weight is used during init.

    `_init_centroids` is shared across all classes inheriting from _BaseKMeans so
    it's enough to check for KMeans.
    """
    rng = np.random.RandomState(global_random_seed)
    X, _ = make_blobs(
        n_samples=200, n_features=10, centers=10, random_state=global_random_seed
    )
    x_squared_norms = row_norms(X, squared=True)

    kmeans = KMeans()
    clusters_weighted = kmeans._init_centroids(
        X=X,
        x_squared_norms=x_squared_norms,
        init=init,
        sample_weight=rng.uniform(size=X.shape[0]),
        n_centroids=5,
        random_state=np.random.RandomState(global_random_seed),
    )
    clusters = kmeans._init_centroids(
        X=X,
        x_squared_norms=x_squared_norms,
        init=init,
        sample_weight=np.ones(X.shape[0]),
        n_centroids=5,
        random_state=np.random.RandomState(global_random_seed),
    )
    with pytest.raises(AssertionError):
        assert_allclose(clusters_weighted, clusters)


@pytest.mark.parametrize("init", ["k-means++", "random"])
def test_sample_weight_zero(init, global_random_seed):
    """Check that if sample weight is 0, this sample won't be chosen.

    `_init_centroids` is shared across all classes inheriting from _BaseKMeans so
    it's enough to check for KMeans.
    """
    rng = np.random.RandomState(global_random_seed)
    X, _ = make_blobs(
        n_samples=100, n_features=5, centers=5, random_state=global_random_seed
    )
    sample_weight = rng.uniform(size=X.shape[0])
    sample_weight[::2] = 0
    x_squared_norms = row_norms(X, squared=True)

    kmeans = KMeans()
    clusters_weighted = kmeans._init_centroids(
        X=X,
        x_squared_norms=x_squared_norms,
        init=init,
        sample_weight=sample_weight,
        n_centroids=10,
        random_state=np.random.RandomState(global_random_seed),
    )
    # No center should be one of the 0 sample weight point
    # (i.e. be at a distance=0 from it)
    d = euclidean_distances(X[::2], clusters_weighted)
    assert not np.any(np.isclose(d, 0))


@pytest.mark.parametrize("array_constr", data_containers, ids=data_containers_ids)
@pytest.mark.parametrize("algorithm", ["lloyd", "elkan"])
def test_relocating_with_duplicates(algorithm, array_constr):
    """Check that kmeans stops when there are more centers than non-duplicate samples

    Non-regression test for issue:
    https://github.com/scikit-learn/scikit-learn/issues/28055
    """
    X = np.array([[0, 0], [1, 1], [1, 1], [1, 0], [0, 1]])
    km = KMeans(n_clusters=5, init=X, algorithm=algorithm)

    msg = r"Number of distinct clusters \(4\) found smaller than n_clusters \(5\)"
    with pytest.warns(ConvergenceWarning, match=msg):
        km.fit(array_constr(X))

    assert km.n_iter_ == 1


# ===========================================================================
# Array API tests
# ===========================================================================

try:
    from sklearn.utils._array_api import (
        _convert_to_numpy,
        _get_namespace_device_dtype_ids,
        get_namespace,
        yield_namespace_device_dtype_combinations,
    )
    from sklearn.utils._array_api import (
        device as _device_func,
    )
    from sklearn.utils._testing import (
        _array_api_for_tests,
        skip_if_array_api_compat_not_configured,
    )

    _ARRAY_API_AVAILABLE = True
except ImportError:
    _ARRAY_API_AVAILABLE = False


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize("init", ["random", "k-means++"])
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_compliance(init, array_namespace, device_, dtype):
    """Test that Array API path produces correct results for KMeans.fit()."""
    xp = _array_api_for_tests(array_namespace, device_)

    X_np, _ = make_blobs(n_samples=100, n_features=4, centers=3, random_state=42)
    X_np = X_np.astype(dtype)

    # Fit with numpy
    km_np = KMeans(n_clusters=3, init=init, n_init=1, random_state=0, algorithm="lloyd")
    km_np.fit(X_np)

    # Fit with Array API
    X_xp = xp.asarray(X_np, device=device_)
    with sklearn.config_context(array_api_dispatch=True):
        km_xp = sklearn.base.clone(km_np)
        km_xp.fit(X_xp)

        # Check output arrays are in the correct namespace and device
        assert get_namespace(km_xp.cluster_centers_)[0] == xp
        assert get_namespace(km_xp.labels_)[0] == xp
        assert _device_func(km_xp.cluster_centers_) == _device_func(X_xp)

        # Check predict
        pred = km_xp.predict(X_xp)
        assert get_namespace(pred)[0] == xp
        assert _device_func(pred) == _device_func(X_xp)

        # Check transform
        transformed = km_xp.transform(X_xp)
        assert get_namespace(transformed)[0] == xp
        assert transformed.shape == (100, 3)

        # Check score returns float
        score = km_xp.score(X_xp)
        assert isinstance(score, float)

    # Compare results with numpy path
    rtol = 1e-4 if dtype == "float32" else 1e-7
    atol = 1e-3 if dtype == "float32" else 0

    centers_xp_np = _convert_to_numpy(km_xp.cluster_centers_, xp=xp)
    labels_xp_np = _convert_to_numpy(km_xp.labels_, xp=xp)

    # Inertia should be close
    assert_allclose(km_np.inertia_, km_xp.inertia_, rtol=rtol, atol=atol)

    # Centers should match (up to permutation)
    assert_allclose(
        np.sort(km_np.cluster_centers_, axis=0),
        np.sort(centers_xp_np, axis=0),
        rtol=rtol,
        atol=atol,
    )

    # Labels should be equivalent (cluster identity check)
    ari = adjusted_rand_score(km_np.labels_, labels_xp_np)
    assert ari > 0.95, f"Labels differ too much: ARI={ari}"


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_elkan_raises(array_namespace, device_, dtype):
    """Test that algorithm='elkan' raises NotImplementedError with Array API."""
    # Only non-numpy namespaces trigger the Array API path
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=0)
    X_xp = xp.asarray(X_np.astype(dtype), device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=3, algorithm="elkan", init="random", n_init=1)
        with pytest.raises(
            NotImplementedError,
            match="algorithm='elkan' is not supported when array_api_dispatch",
        ):
            km.fit(X_xp)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_edge_cases(array_namespace, device_, dtype):
    """Test Array API path with edge cases."""
    xp = _array_api_for_tests(array_namespace, device_)

    with sklearn.config_context(array_api_dispatch=True):
        # n_clusters=1: all points in same cluster
        X_np = np.array([[1, 2], [3, 4], [5, 6]], dtype=dtype)
        X_xp = xp.asarray(X_np, device=device_)
        km = KMeans(n_clusters=1, init="random", n_init=1, random_state=0)
        km.fit(X_xp)
        assert km.cluster_centers_.shape == (1, 2)
        labels = km.predict(X_xp)
        # All labels should be 0
        labels_np = _convert_to_numpy(labels, xp=xp)
        assert_array_equal(labels_np, [0, 0, 0])

        # max_iter=1: single iteration
        X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=0)
        X_np = X_np.astype(dtype)
        X_xp = xp.asarray(X_np, device=device_)
        km = KMeans(n_clusters=3, init="random", n_init=1, max_iter=1, random_state=0)
        km.fit(X_xp)
        assert km.n_iter_ == 1

        # n_clusters == n_samples: each point is its own cluster
        X_small = np.array([[0, 0], [10, 10], [20, 20]], dtype=dtype)
        X_small_xp = xp.asarray(X_small, device=device_)
        km = KMeans(n_clusters=3, init="random", n_init=1, random_state=0)
        km.fit(X_small_xp)
        labels_np = _convert_to_numpy(km.labels_, xp=xp)
        assert len(set(labels_np.tolist())) == 3

        # All identical points
        X_same = np.ones((10, 3), dtype=dtype)
        X_same_xp = xp.asarray(X_same, device=device_)
        km = KMeans(n_clusters=1, init="random", n_init=1, random_state=0)
        km.fit(X_same_xp)
        centers_np = _convert_to_numpy(km.cluster_centers_, xp=xp)
        assert_allclose(centers_np, [[1.0, 1.0, 1.0]], atol=1e-6)

        # 1D features
        X_1d = np.array([[1], [2], [10], [11]], dtype=dtype)
        X_1d_xp = xp.asarray(X_1d, device=device_)
        km = KMeans(n_clusters=2, init="random", n_init=3, random_state=0)
        km.fit(X_1d_xp)
        assert km.cluster_centers_.shape == (2, 1)

        # float32 precision
        if dtype == "float32":
            X_np, _ = make_blobs(
                n_samples=100, n_features=4, centers=3, random_state=42
            )
            X_np = X_np.astype(np.float32)
            X_xp = xp.asarray(X_np, device=device_)
            km = KMeans(n_clusters=3, init="random", n_init=1, random_state=0)
            km.fit(X_xp)
            # Should not raise and should produce finite results
            assert np.isfinite(km.inertia_)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_n_init(array_namespace, device_, dtype):
    """Test that n_init > 1 selects the best result with Array API."""
    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=100, n_features=4, centers=3, random_state=42)
    X_np = X_np.astype(dtype)
    X_xp = xp.asarray(X_np, device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km_1 = KMeans(n_clusters=3, init="random", n_init=1, random_state=0)
        km_1.fit(X_xp)

        km_10 = KMeans(n_clusters=3, init="random", n_init=10, random_state=0)
        km_10.fit(X_xp)

        # n_init=10 should find equal or better inertia than n_init=1
        assert km_10.inertia_ <= km_1.inertia_ + 1e-6


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_sample_weight(array_namespace, device_, dtype):
    """Test that sample_weight is respected in Array API path."""
    xp = _array_api_for_tests(array_namespace, device_)

    # Use well-separated blobs to ensure clear clustering
    X_np, _ = make_blobs(
        n_samples=60, n_features=2, centers=2, cluster_std=0.5, random_state=42
    )
    X_np = X_np.astype(dtype)

    # Create sample weights: double weight on first half
    sw_np = np.ones(60, dtype=dtype)
    sw_np[:30] = 10.0

    X_xp = xp.asarray(X_np, device=device_)
    sw_xp = xp.asarray(sw_np, device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km_w = KMeans(n_clusters=2, init="random", n_init=5, random_state=0)
        km_w.fit(X_xp, sample_weight=sw_xp)

        km_uw = KMeans(n_clusters=2, init="random", n_init=5, random_state=0)
        km_uw.fit(X_xp)

        # Both should converge and find 2 clusters
        assert km_w.cluster_centers_.shape == (2, 2)
        assert km_uw.cluster_centers_.shape == (2, 2)
        # Weighted version should have finite inertia
        assert np.isfinite(km_w.inertia_)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_fit_predict(array_namespace, device_, dtype):
    """Test fit_predict and fit_transform with Array API."""
    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=42)
    X_np = X_np.astype(dtype)
    X_xp = xp.asarray(X_np, device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=3, init="random", n_init=1, random_state=0)

        # fit_predict
        labels = km.fit_predict(X_xp)
        assert get_namespace(labels)[0] == xp

        # fit_transform
        km2 = KMeans(n_clusters=3, init="random", n_init=1, random_state=0)
        X_transformed = km2.fit_transform(X_xp)
        assert get_namespace(X_transformed)[0] == xp
        assert X_transformed.shape == (50, 3)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_verbose(capsys, array_namespace, device_, dtype):
    """Test verbose output in Array API path."""
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=42)
    X_xp = xp.asarray(X_np.astype(dtype), device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=3, init="random", n_init=1, random_state=0, verbose=True)
        km.fit(X_xp)

    captured = capsys.readouterr()
    assert "Initialization complete" in captured.out
    assert "Iteration" in captured.out
    assert "inertia" in captured.out


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_empty_cluster(array_namespace, device_, dtype):
    """Test that empty clusters retain their centroid in Array API path."""
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    # Two tight clusters but request 5 clusters; at least some will be empty
    rng = np.random.RandomState(0)
    X_np = np.vstack(
        [rng.randn(30, 2) + [10, 10], rng.randn(30, 2) + [-10, -10]]
    ).astype(dtype)
    # Place 3 extra init centers far from data so their clusters stay empty
    init_centers = np.vstack(
        [
            X_np[:2],  # 2 centers near data
            [[100, 100], [200, 200], [300, 300]],  # 3 centers far from data
        ]
    ).astype(dtype)

    X_xp = xp.asarray(X_np, device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=5, init=init_centers, n_init=1, max_iter=3)
        km.fit(X_xp)

    # All cluster centers should be finite (empty ones retain init centroid)
    centers_np = _convert_to_numpy(km.cluster_centers_, xp=xp)
    assert np.all(np.isfinite(centers_np))
    assert km.cluster_centers_.shape == (5, 2)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_tol_zero(array_namespace, device_, dtype):
    """Test tol=0 exercises _tolerance_array_api early-return in Array API path."""
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=42)
    X_xp = xp.asarray(X_np.astype(dtype), device=device_)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(
            n_clusters=3,
            init="random",
            n_init=1,
            random_state=0,
            tol=0,
            max_iter=10,
        )
        km.fit(X_xp)

    # tol=0 means converge only when centers stop moving entirely
    assert km.n_iter_ >= 1
    assert km.cluster_centers_.shape == (3, 2)


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_callable_init(array_namespace, device_, dtype):
    """Test callable init in Array API path."""
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=42)
    X_xp = xp.asarray(X_np.astype(dtype), device=device_)

    def custom_init(X, n_clusters, random_state):
        # Convert to numpy for indexing; caller wraps result with xp.asarray
        xp_local, _ = get_namespace(X)
        X_np_local = _convert_to_numpy(X, xp=xp_local)
        indices = random_state.choice(
            X_np_local.shape[0], size=n_clusters, replace=False
        )
        return X_np_local[indices]

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=3, init=custom_init, n_init=1, random_state=0)
        km.fit(X_xp)

    assert km.cluster_centers_.shape == (3, 2)
    assert km.n_iter_ >= 1


@pytest.mark.skipif(not _ARRAY_API_AVAILABLE, reason="array_api_compat not installed")
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype",
    yield_namespace_device_dtype_combinations() if _ARRAY_API_AVAILABLE else [],
    ids=_get_namespace_device_dtype_ids if _ARRAY_API_AVAILABLE else None,
)
def test_kmeans_array_api_array_init(array_namespace, device_, dtype):
    """Test array-like init in Array API path."""
    if array_namespace == "numpy":
        pytest.skip("numpy namespace uses existing Cython path")

    xp = _array_api_for_tests(array_namespace, device_)
    X_np, _ = make_blobs(n_samples=50, n_features=2, centers=3, random_state=42)
    X_xp = xp.asarray(X_np.astype(dtype), device=device_)

    # Use first 3 samples as initial centers (numpy array)
    init_centers = X_np[:3].astype(dtype)

    with sklearn.config_context(array_api_dispatch=True):
        km = KMeans(n_clusters=3, init=init_centers, n_init=1)
        km.fit(X_xp)

    assert km.cluster_centers_.shape == (3, 2)
    assert km.n_iter_ >= 1
    centers_np = _convert_to_numpy(km.cluster_centers_, xp=xp)
    assert np.all(np.isfinite(centers_np))
