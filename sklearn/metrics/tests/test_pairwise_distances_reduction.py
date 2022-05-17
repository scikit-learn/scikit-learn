import numpy as np
import pytest
import threadpoolctl
from numpy.testing import assert_array_equal, assert_allclose
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist

from sklearn.metrics._pairwise_distances_reduction import (
    PairwiseDistancesReduction,
    PairwiseDistancesArgKmin,
    PairwiseDistancesRadiusNeighborhood,
    _sqeuclidean_row_norms,
)

from sklearn.metrics import euclidean_distances
from sklearn.utils.fixes import sp_version, parse_version

# Common supported metric between scipy.spatial.distance.cdist
# and PairwiseDistancesReduction.
# This allows constructing tests to check consistency of results
# of concrete PairwiseDistancesReduction on some metrics using APIs
# from scipy and numpy.
CDIST_PAIRWISE_DISTANCES_REDUCTION_COMMON_METRICS = [
    "braycurtis",
    "canberra",
    "chebyshev",
    "cityblock",
    "euclidean",
    "minkowski",
    "seuclidean",
]


def _get_metric_params_list(metric: str, n_features: int, seed: int = 1):
    """Return list of dummy DistanceMetric kwargs for tests."""

    # Distinguishing on cases not to compute unneeded datastructures.
    rng = np.random.RandomState(seed)

    if metric == "minkowski":
        minkowski_kwargs = [dict(p=1.5), dict(p=2), dict(p=3), dict(p=np.inf)]
        if sp_version >= parse_version("1.8.0.dev0"):
            # TODO: remove the test once we no longer support scipy < 1.8.0.
            # Recent scipy versions accept weights in the Minkowski metric directly:
            # type: ignore
            minkowski_kwargs.append(dict(p=3, w=rng.rand(n_features)))

        return minkowski_kwargs

    # TODO: remove this case for "wminkowski" once we no longer support scipy < 1.8.0.
    if metric == "wminkowski":
        weights = rng.random_sample(n_features)
        weights /= weights.sum()
        wminkowski_kwargs = [dict(p=1.5, w=weights)]
        if sp_version < parse_version("1.8.0.dev0"):
            # wminkowski was removed in scipy 1.8.0 but should work for previous
            # versions.
            wminkowski_kwargs.append(dict(p=3, w=rng.rand(n_features)))
        return wminkowski_kwargs

    if metric == "seuclidean":
        return [dict(V=rng.rand(n_features))]

    # Case of: "euclidean", "manhattan", "chebyshev", "haversine" or any other metric.
    # In those cases, no kwargs is needed.
    return [{}]


def assert_argkmin_results_equality(ref_dist, dist, ref_indices, indices):
    assert_array_equal(
        ref_indices,
        indices,
        err_msg="Query vectors have different neighbors' indices",
    )
    assert_allclose(
        ref_dist,
        dist,
        err_msg="Query vectors have different neighbors' distances",
        rtol=1e-7,
    )


def assert_radius_neighborhood_results_equality(ref_dist, dist, ref_indices, indices):
    # We get arrays of arrays and we need to check for individual pairs
    for i in range(ref_dist.shape[0]):
        assert_array_equal(
            ref_indices[i],
            indices[i],
            err_msg=f"Query vector #{i} has different neighbors' indices",
        )
        assert_allclose(
            ref_dist[i],
            dist[i],
            err_msg=f"Query vector #{i} has different neighbors' distances",
            rtol=1e-7,
        )


ASSERT_RESULT = {
    PairwiseDistancesArgKmin: assert_argkmin_results_equality,
    PairwiseDistancesRadiusNeighborhood: assert_radius_neighborhood_results_equality,
}


def test_pairwise_distances_reduction_is_usable_for():
    rng = np.random.RandomState(0)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    metric = "euclidean"
    assert PairwiseDistancesReduction.is_usable_for(X, Y, metric)
    assert not PairwiseDistancesReduction.is_usable_for(
        X.astype(np.int64), Y.astype(np.int64), metric
    )

    assert not PairwiseDistancesReduction.is_usable_for(X, Y, metric="pyfunc")
    # TODO: remove once 32 bits datasets are supported
    assert not PairwiseDistancesReduction.is_usable_for(X.astype(np.float32), Y, metric)
    assert not PairwiseDistancesReduction.is_usable_for(X, Y.astype(np.int32), metric)

    # TODO: remove once sparse matrices are supported
    assert not PairwiseDistancesReduction.is_usable_for(csr_matrix(X), Y, metric)
    assert not PairwiseDistancesReduction.is_usable_for(X, csr_matrix(Y), metric)


def test_argkmin_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    k = 5
    metric = "euclidean"

    msg = (
        "Only 64bit float datasets are supported at this time, "
        "got: X.dtype=float32 and Y.dtype=float64"
    )
    with pytest.raises(ValueError, match=msg):
        PairwiseDistancesArgKmin.compute(
            X=X.astype(np.float32), Y=Y, k=k, metric=metric
        )

    msg = (
        "Only 64bit float datasets are supported at this time, "
        "got: X.dtype=float64 and Y.dtype=int32"
    )
    with pytest.raises(ValueError, match=msg):
        PairwiseDistancesArgKmin.compute(X=X, Y=Y.astype(np.int32), k=k, metric=metric)

    with pytest.raises(ValueError, match="k == -1, must be >= 1."):
        PairwiseDistancesArgKmin.compute(X=X, Y=Y, k=-1, metric=metric)

    with pytest.raises(ValueError, match="k == 0, must be >= 1."):
        PairwiseDistancesArgKmin.compute(X=X, Y=Y, k=0, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        PairwiseDistancesArgKmin.compute(X=X, Y=Y, k=k, metric="wrong metric")

    with pytest.raises(
        ValueError, match=r"Buffer has wrong number of dimensions \(expected 2, got 1\)"
    ):
        PairwiseDistancesArgKmin.compute(
            X=np.array([1.0, 2.0]), Y=Y, k=k, metric=metric
        )

    with pytest.raises(ValueError, match="ndarray is not C-contiguous"):
        PairwiseDistancesArgKmin.compute(
            X=np.asfortranarray(X), Y=Y, k=k, metric=metric
        )

    unused_metric_kwargs = {"p": 3}

    message = (
        r"Some metric_kwargs have been passed \({'p': 3}\) but aren't usable for this"
        r" case \("
        r"FastEuclideanPairwiseDistancesArgKmin\) and will be ignored."
    )

    with pytest.warns(UserWarning, match=message):
        PairwiseDistancesArgKmin.compute(
            X=X, Y=Y, k=k, metric=metric, metric_kwargs=unused_metric_kwargs
        )


def test_radius_neighborhood_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    radius = 5
    metric = "euclidean"

    with pytest.raises(
        ValueError,
        match=(
            "Only 64bit float datasets are supported at this time, "
            "got: X.dtype=float32 and Y.dtype=float64"
        ),
    ):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=X.astype(np.float32), Y=Y, radius=radius, metric=metric
        )

    with pytest.raises(
        ValueError,
        match=(
            "Only 64bit float datasets are supported at this time, "
            "got: X.dtype=float64 and Y.dtype=int32"
        ),
    ):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=X, Y=Y.astype(np.int32), radius=radius, metric=metric
        )

    with pytest.raises(ValueError, match="radius == -1.0, must be >= 0."):
        PairwiseDistancesRadiusNeighborhood.compute(X=X, Y=Y, radius=-1, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=X, Y=Y, radius=radius, metric="wrong metric"
        )

    with pytest.raises(
        ValueError, match=r"Buffer has wrong number of dimensions \(expected 2, got 1\)"
    ):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=np.array([1.0, 2.0]), Y=Y, radius=radius, metric=metric
        )

    with pytest.raises(ValueError, match="ndarray is not C-contiguous"):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=np.asfortranarray(X), Y=Y, radius=radius, metric=metric
        )

    unused_metric_kwargs = {"p": 3}

    message = (
        r"Some metric_kwargs have been passed \({'p': 3}\) but aren't usable for this"
        r" case \(FastEuclideanPairwiseDistancesRadiusNeighborhood\) and will be"
        r" ignored."
    )

    with pytest.warns(UserWarning, match=message):
        PairwiseDistancesRadiusNeighborhood.compute(
            X=X, Y=Y, radius=radius, metric=metric, metric_kwargs=unused_metric_kwargs
        )


@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin, PairwiseDistancesRadiusNeighborhood],
)
def test_chunk_size_agnosticism(
    global_random_seed,
    PairwiseDistancesReduction,
    n_samples,
    chunk_size,
    n_features=100,
    dtype=np.float64,
):
    # Results should not depend on the chunk size
    rng = np.random.RandomState(global_random_seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    parameter = (
        10
        if PairwiseDistancesReduction is PairwiseDistancesArgKmin
        # Scaling the radius slightly with the numbers of dimensions
        else 10 ** np.log(n_features)
    )

    ref_dist, ref_indices = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        return_distance=True,
    )

    dist, indices = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        chunk_size=chunk_size,
        return_distance=True,
    )

    ASSERT_RESULT[PairwiseDistancesReduction](ref_dist, dist, ref_indices, indices)


@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin, PairwiseDistancesRadiusNeighborhood],
)
def test_n_threads_agnosticism(
    global_random_seed,
    PairwiseDistancesReduction,
    n_samples,
    chunk_size,
    n_features=100,
    dtype=np.float64,
):
    # Results should not depend on the number of threads
    rng = np.random.RandomState(global_random_seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    parameter = (
        10
        if PairwiseDistancesReduction is PairwiseDistancesArgKmin
        # Scaling the radius slightly with the numbers of dimensions
        else 10 ** np.log(n_features)
    )

    ref_dist, ref_indices = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        return_distance=True,
    )

    with threadpoolctl.threadpool_limits(limits=1, user_api="openmp"):
        dist, indices = PairwiseDistancesReduction.compute(
            X, Y, parameter, return_distance=True
        )

    ASSERT_RESULT[PairwiseDistancesReduction](ref_dist, dist, ref_indices, indices)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("metric", PairwiseDistancesReduction.valid_metrics())
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin, PairwiseDistancesRadiusNeighborhood],
)
def test_strategies_consistency(
    global_random_seed,
    PairwiseDistancesReduction,
    metric,
    n_samples,
    n_features=10,
    dtype=np.float64,
):

    rng = np.random.RandomState(global_random_seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    # Haversine distance only accepts 2D data
    if metric == "haversine":
        X = np.ascontiguousarray(X[:, :2])
        Y = np.ascontiguousarray(Y[:, :2])

    parameter = (
        10
        if PairwiseDistancesReduction is PairwiseDistancesArgKmin
        # Scaling the radius slightly with the numbers of dimensions
        else 10 ** np.log(n_features)
    )

    dist_par_X, indices_par_X = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        metric=metric,
        # Taking the first
        metric_kwargs=_get_metric_params_list(
            metric, n_features, seed=global_random_seed
        )[0],
        # To be sure to use parallelization
        chunk_size=n_samples // 4,
        strategy="parallel_on_X",
        return_distance=True,
    )

    dist_par_Y, indices_par_Y = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        metric=metric,
        # Taking the first
        metric_kwargs=_get_metric_params_list(
            metric, n_features, seed=global_random_seed
        )[0],
        # To be sure to use parallelization
        chunk_size=n_samples // 4,
        strategy="parallel_on_Y",
        return_distance=True,
    )

    ASSERT_RESULT[PairwiseDistancesReduction](
        dist_par_X,
        dist_par_Y,
        indices_par_X,
        indices_par_Y,
    )


# "Concrete PairwiseDistancesReductions"-specific tests

# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("n_features", [50, 500])
@pytest.mark.parametrize("translation", [0, 1e6])
@pytest.mark.parametrize("metric", CDIST_PAIRWISE_DISTANCES_REDUCTION_COMMON_METRICS)
@pytest.mark.parametrize("strategy", ("parallel_on_X", "parallel_on_Y"))
def test_pairwise_distances_argkmin(
    global_random_seed,
    n_features,
    translation,
    metric,
    strategy,
    n_samples=100,
    k=10,
    dtype=np.float64,
):
    rng = np.random.RandomState(global_random_seed)
    spread = 1000
    X = translation + rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = translation + rng.rand(n_samples, n_features).astype(dtype) * spread

    # Haversine distance only accepts 2D data
    if metric == "haversine":
        X = np.ascontiguousarray(X[:, :2])
        Y = np.ascontiguousarray(Y[:, :2])

    metric_kwargs = _get_metric_params_list(metric, n_features)[0]

    # Reference for argkmin results
    if metric == "euclidean":
        # Compare to scikit-learn GEMM optimized implementation
        dist_matrix = euclidean_distances(X, Y)
    else:
        dist_matrix = cdist(X, Y, metric=metric, **metric_kwargs)
    # Taking argkmin (indices of the k smallest values)
    argkmin_indices_ref = np.argsort(dist_matrix, axis=1)[:, :k]
    # Getting the associated distances
    argkmin_distances_ref = np.zeros(argkmin_indices_ref.shape, dtype=np.float64)
    for row_idx in range(argkmin_indices_ref.shape[0]):
        argkmin_distances_ref[row_idx] = dist_matrix[
            row_idx, argkmin_indices_ref[row_idx]
        ]

    argkmin_distances, argkmin_indices = PairwiseDistancesArgKmin.compute(
        X,
        Y,
        k,
        metric=metric,
        metric_kwargs=metric_kwargs,
        return_distance=True,
        # So as to have more than a chunk, forcing parallelism.
        chunk_size=n_samples // 4,
        strategy=strategy,
    )

    ASSERT_RESULT[PairwiseDistancesArgKmin](
        argkmin_distances, argkmin_distances_ref, argkmin_indices, argkmin_indices_ref
    )


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("n_features", [50, 500])
@pytest.mark.parametrize("translation", [0, 1e6])
@pytest.mark.parametrize("metric", CDIST_PAIRWISE_DISTANCES_REDUCTION_COMMON_METRICS)
@pytest.mark.parametrize("strategy", ("parallel_on_X", "parallel_on_Y"))
def test_pairwise_distances_radius_neighbors(
    global_random_seed,
    n_features,
    translation,
    metric,
    strategy,
    n_samples=100,
    dtype=np.float64,
):
    rng = np.random.RandomState(global_random_seed)
    spread = 1000
    radius = spread * np.log(n_features)
    X = translation + rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = translation + rng.rand(n_samples, n_features).astype(dtype) * spread

    metric_kwargs = _get_metric_params_list(
        metric, n_features, seed=global_random_seed
    )[0]

    # Reference for argkmin results
    if metric == "euclidean":
        # Compare to scikit-learn GEMM optimized implementation
        dist_matrix = euclidean_distances(X, Y)
    else:
        dist_matrix = cdist(X, Y, metric=metric, **metric_kwargs)

    # Getting the neighbors for a given radius
    neigh_indices_ref = []
    neigh_distances_ref = []

    for row in dist_matrix:
        ind = np.arange(row.shape[0])[row <= radius]
        dist = row[ind]

        sort = np.argsort(dist)
        ind, dist = ind[sort], dist[sort]

        neigh_indices_ref.append(ind)
        neigh_distances_ref.append(dist)

    neigh_indices_ref = np.array(neigh_indices_ref)
    neigh_distances_ref = np.array(neigh_distances_ref)

    neigh_distances, neigh_indices = PairwiseDistancesRadiusNeighborhood.compute(
        X,
        Y,
        radius,
        metric=metric,
        metric_kwargs=metric_kwargs,
        return_distance=True,
        # So as to have more than a chunk, forcing parallelism.
        chunk_size=n_samples // 4,
        strategy=strategy,
        sort_results=True,
    )

    ASSERT_RESULT[PairwiseDistancesRadiusNeighborhood](
        neigh_distances, neigh_distances_ref, neigh_indices, neigh_indices_ref
    )


@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("n_features", [5, 10, 100])
@pytest.mark.parametrize("num_threads", [1, 2, 8])
def test_sqeuclidean_row_norms(
    global_random_seed,
    n_samples,
    n_features,
    num_threads,
    dtype=np.float64,
):
    rng = np.random.RandomState(global_random_seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread

    sq_row_norm_reference = np.linalg.norm(X, axis=1) ** 2
    sq_row_norm = np.asarray(_sqeuclidean_row_norms(X, num_threads=num_threads))

    assert_allclose(sq_row_norm_reference, sq_row_norm)
