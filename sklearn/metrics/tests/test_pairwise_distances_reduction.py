import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_allclose
from scipy.sparse import csr_matrix

from sklearn.metrics._pairwise_distances_reduction import (
    PairwiseDistancesReduction,
    PairwiseDistancesArgKmin,
    _sqeuclidean_row_norms,
)

from sklearn.utils import _in_unstable_openblas_configuration
from sklearn.utils.fixes import sp_version, parse_version
from sklearn.utils._testing import fails_if_unstable_openblas


def _get_dummy_metric_params_list(metric: str, n_features: int):
    """Return list of dummy DistanceMetric kwargs for tests."""

    # Distinguishing on cases not to compute unneeded datastructures.
    rng = np.random.RandomState(1)

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

    if metric == "mahalanobis":
        V = rng.random_sample((n_features, n_features))
        # This makes VI is positive-semidefinite, which is a
        # necessary condition to get nonsingular precision matrix.
        VI = np.dot(V, V.T) + 3 * np.eye(n_features)

        return [dict(VI=VI)]

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


ASSERT_RESULT = {
    PairwiseDistancesArgKmin: assert_argkmin_results_equality,
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

    assert not PairwiseDistancesReduction.is_usable_for(X[0], Y, metric)
    assert not PairwiseDistancesReduction.is_usable_for(X, Y[0], metric)

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

    with pytest.raises(
        ValueError, match="Only 64bit float datasets are supported for X and Y."
    ):
        PairwiseDistancesArgKmin.compute(
            X=X.astype(np.float32), Y=Y, k=k, metric=metric
        )

    with pytest.raises(
        ValueError, match="Only 64bit float datasets are supported for X and Y."
    ):
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


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin],
)
def test_chunk_size_agnosticism(
    PairwiseDistancesReduction,
    seed,
    n_samples,
    chunk_size,
    n_features=100,
    dtype=np.float64,
):
    # Results should not depend on the chunk size
    rng = np.random.RandomState(seed)
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


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin],
)
def test_n_threads_agnosticism(
    PairwiseDistancesReduction,
    seed,
    n_samples,
    chunk_size,
    n_features=100,
    dtype=np.float64,
):
    # Results should not depend on the number of threads
    rng = np.random.RandomState(seed)
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
        X, Y, parameter, n_threads=1, return_distance=True
    )

    ASSERT_RESULT[PairwiseDistancesReduction](ref_dist, dist, ref_indices, indices)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("metric", PairwiseDistancesReduction.valid_metrics())
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin],
)
def test_strategies_consistency(
    PairwiseDistancesReduction,
    metric,
    n_samples,
    seed,
    n_features=10,
    dtype=np.float64,
):
    # Results obtained using both parallelization strategies must be identical
    if _in_unstable_openblas_configuration() and metric in ("sqeuclidean", "euclidean"):
        pytest.xfail(
            "OpenBLAS (used for '(sq)euclidean') is unstable in this configuration"
        )

    rng = np.random.RandomState(seed)
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
        metric_kwargs=_get_dummy_metric_params_list(metric, n_features)[0],
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
        metric_kwargs=_get_dummy_metric_params_list(metric, n_features)[0],
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


@fails_if_unstable_openblas
@pytest.mark.parametrize("n_features", [50, 500])
@pytest.mark.parametrize("translation", [10 ** i for i in [4, 8]])
@pytest.mark.parametrize("metric", PairwiseDistancesReduction.valid_metrics())
@pytest.mark.parametrize(
    "PairwiseDistancesReduction",
    [PairwiseDistancesArgKmin],
)
def test_euclidean_translation_invariance(
    n_features,
    translation,
    metric,
    PairwiseDistancesReduction,
    n_samples=1000,
    dtype=np.float64,
):
    # The reduction must be translation invariant.
    parameter = (
        10
        if PairwiseDistancesReduction is PairwiseDistancesArgKmin
        # Scaling the radius slightly with the numbers of dimensions
        else 10 ** np.log(n_features)
    )

    rng = np.random.RandomState(0)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    # Haversine distance only accepts 2D data
    if metric == "haversine":
        X = np.ascontiguousarray(X[:, :2])
        Y = np.ascontiguousarray(Y[:, :2])

    reference_dist, reference_indices = PairwiseDistancesReduction.compute(
        X,
        Y,
        parameter,
        metric=metric,
        metric_kwargs=_get_dummy_metric_params_list(metric, n_features)[0],
        return_distance=True,
    )

    dist, indices = PairwiseDistancesReduction.compute(
        X + 0,
        Y + 0,
        parameter,
        metric=metric,
        metric_kwargs=_get_dummy_metric_params_list(metric, n_features)[0],
        return_distance=True,
    )

    ASSERT_RESULT[PairwiseDistancesReduction](
        reference_dist, dist, reference_indices, indices
    )


@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("n_samples", [100, 1000])
@pytest.mark.parametrize("n_features", [5, 10, 100])
@pytest.mark.parametrize("num_threads", [1, 2, 8])
def test_sqeuclidean_row_norms(
    seed,
    n_samples,
    n_features,
    num_threads,
    dtype=np.float64,
):
    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread

    sq_row_norm_reference = np.linalg.norm(X, axis=1) ** 2
    sq_row_norm = np.asarray(_sqeuclidean_row_norms(X, num_threads=num_threads))

    assert_allclose(sq_row_norm_reference, sq_row_norm)
