import itertools

import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_allclose
from scipy.sparse import (
    bsr_matrix,
    coo_matrix,
    csc_matrix,
    csr_matrix,
    dia_matrix,
    dok_matrix,
    lil_matrix,
)

from sklearn.metrics._dist_metrics import (
    DenseDenseDatasetsPair,
    DenseSparseDatasetsPair,
    SparseDenseDatasetsPair,
    SparseSparseDatasetsPair,
)

from sklearn.metrics._pairwise_distances_reduction import (
    PairwiseDistancesReduction,
    ArgKmin,
    RadiusNeighborhood,
    FastSquaredEuclideanArgKmin,
    FastSquaredEuclideanRadiusNeighborhood,
)

from sklearn.utils import _in_unstable_openblas_configuration
from sklearn.utils._testing import (
    fails_if_unstable_openblas,
    get_dummy_metric_kwargs,
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


def assert_argkmin_results_equality(ref_dist, dist, ref_indices, indices):
    assert_array_equal(
        ref_indices,
        indices,
        err_msg="Query vectors have different neighbors' indices",
    )
    assert_allclose(
        ref_dist,
        dist,
        err_msg="Query vectors havehas different neighbors' distances",
        rtol=1e-7,
    )


def test_pairwise_distances_reduction_is_usable_for():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    metric = "euclidean"
    assert PairwiseDistancesReduction.is_usable_for(X, Y, metric)
    assert PairwiseDistancesReduction.is_usable_for(
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
    assert not PairwiseDistancesReduction.is_usable_for(X, csc_matrix(Y), metric)


def test_argkmin_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    k = 5
    metric = "euclidean"

    with pytest.raises(
        ValueError, match="32bits datasets aren't supported for X and Y yet."
    ):
        ArgKmin.get_for(X=X.astype(np.float32), Y=Y, k=k, metric=metric)

    with pytest.raises(
        ValueError, match="32bits datasets aren't supported for X and Y yet."
    ):
        ArgKmin.get_for(X=X, Y=Y.astype(np.int32), k=k, metric=metric)

    with pytest.raises(ValueError, match="k == -1, must be >= 1."):
        ArgKmin.get_for(X=X, Y=Y, k=-1, metric=metric)

    with pytest.raises(ValueError, match="k == 0, must be >= 1."):
        ArgKmin.get_for(X=X, Y=Y, k=0.1, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        ArgKmin.get_for(X=X, Y=Y, k=k, metric="wrong metric")

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        ArgKmin.get_for(X=np.array([1.0, 2.0]), Y=Y, k=k, metric=metric)

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        ArgKmin.get_for(X=X, Y=np.array([1.0, 2.0]), k=k, metric=metric)

    with pytest.raises(
        ValueError, match="Vectors of X and Y must have the same number of dimensions"
    ):
        ArgKmin.get_for(X=X[:, ::2], Y=Y, k=k, metric=metric)


def test_radius_neighborhood_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    radius = 5
    metric = "euclidean"

    with pytest.raises(
        ValueError, match="32bits datasets aren't supported for X and Y yet."
    ):
        RadiusNeighborhood.get_for(
            X=X.astype(np.float32), Y=Y, radius=radius, metric=metric
        )

    with pytest.raises(
        ValueError, match="32bits datasets aren't supported for X and Y yet."
    ):
        RadiusNeighborhood.get_for(
            X=X, Y=Y.astype(np.int32), radius=radius, metric=metric
        )

    with pytest.raises(ValueError, match="radius == -1.0, must be >= 0."):
        RadiusNeighborhood.get_for(X=X, Y=Y, radius=-1, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        RadiusNeighborhood.get_for(X=X, Y=Y, radius=radius, metric="wrong metric")

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        RadiusNeighborhood.get_for(
            X=np.array([1.0, 2.0]), Y=Y, radius=radius, metric=metric
        )

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        RadiusNeighborhood.get_for(
            X=X, Y=np.array([1.0, 2.0]), radius=radius, metric=metric
        )

    with pytest.raises(
        ValueError, match="Vectors of X and Y must have the same number of dimensions"
    ):
        RadiusNeighborhood.get_for(X=X[:, ::2], Y=Y, radius=radius, metric=metric)


@fails_if_unstable_openblas
@pytest.mark.filterwarnings("ignore:Constructing a DIA matrix")
@pytest.mark.parametrize(
    "PairwiseDistancesReduction, FastSquaredPairwiseDistancesReduction",
    [
        (ArgKmin, FastSquaredEuclideanArgKmin),
        (RadiusNeighborhood, FastSquaredEuclideanRadiusNeighborhood),
    ],
)
def test_pairwise_distances_reduction_factory_method(
    PairwiseDistancesReduction, FastSquaredPairwiseDistancesReduction
):
    # Test all the combinations of DatasetsPair for creation
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    metric = "euclidean"

    # Dummy value for k or radius
    dummy_arg = 5

    dense_dense_instance = PairwiseDistancesReduction.get_for(X, Y, dummy_arg, metric)
    assert isinstance(dense_dense_instance.datasets_pair, DenseDenseDatasetsPair)

    sparse_matrix_constructors = [
        lil_matrix,
        csc_matrix,
        csr_matrix,
        bsr_matrix,
        coo_matrix,
        dia_matrix,
        dok_matrix,
    ]

    for c_X, c_Y in itertools.combinations_with_replacement(
        sparse_matrix_constructors, r=2
    ):
        sparse_sparse_instance = PairwiseDistancesReduction.get_for(
            c_X(X), c_Y(Y), dummy_arg, metric
        )
        assert isinstance(
            sparse_sparse_instance.datasets_pair, SparseSparseDatasetsPair
        )

    for constructor in sparse_matrix_constructors:
        dense_sparse_instance = PairwiseDistancesReduction.get_for(
            X, constructor(Y), dummy_arg, metric=metric
        )
        assert isinstance(dense_sparse_instance.datasets_pair, DenseSparseDatasetsPair)

        sparse_dense_instance = PairwiseDistancesReduction.get_for(
            constructor(X), Y, dummy_arg, metric=metric
        )
        assert isinstance(sparse_dense_instance.datasets_pair, SparseDenseDatasetsPair)

    # Test specialisations creation
    fast_sqeuclidean_instance = PairwiseDistancesReduction.get_for(
        X, Y, dummy_arg, metric="fast_sqeuclidean"
    )
    assert isinstance(fast_sqeuclidean_instance, PairwiseDistancesReduction)
    assert isinstance(fast_sqeuclidean_instance, FastSquaredPairwiseDistancesReduction)


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
def test_argkmin_chunk_size_agnosticism(
    seed,
    n_samples,
    chunk_size,
    k=10,
    metric="fast_sqeuclidean",
    n_features=100,
    dtype=np.float64,
):
    # ArgKmin results should not depend on the chunk size
    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    ref_dist, ref_indices = ArgKmin.get_for(X, Y, k=k, metric="euclidean").compute(
        return_distance=True
    )

    dist, indices = ArgKmin.get_for(
        X, Y, k=k, metric=metric, chunk_size=chunk_size
    ).compute(return_distance=True)

    assert_argkmin_results_equality(ref_dist, dist, ref_indices, indices)


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("chunk_size", [50, 512, 1024])
def test_radius_neighborhood_chunk_size_agnosticism(
    seed,
    n_samples,
    chunk_size,
    radius=10.0,
    metric="fast_sqeuclidean",
    n_features=100,
    dtype=np.float64,
):
    # RadiusNeighborhood results should not depend on the chunk size
    rng = np.random.RandomState(seed)
    spread = 100

    # Scaling the radius with the dimensions
    scaled_radius = radius * np.log(n_features)
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    ref_dist, ref_indices = RadiusNeighborhood.get_for(
        X, Y, radius=scaled_radius, metric="euclidean"
    ).compute(return_distance=True)

    dist, indices = RadiusNeighborhood.get_for(
        X, Y, radius=scaled_radius, metric=metric, chunk_size=chunk_size
    ).compute(return_distance=True)

    assert_radius_neighborhood_results_equality(ref_dist, dist, ref_indices, indices)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("metric", RadiusNeighborhood.valid_metrics())
def test_argkmin_strategies_consistency(
    metric,
    n_samples,
    seed,
    n_features=10,
    k=10,
    dtype=np.float64,
):
    # ArgKmin results obtained using both parallelization strategies
    # must be identical
    if _in_unstable_openblas_configuration() and metric == "fast_sqeuclidean":
        pytest.xfail(
            "OpenBLAS (used for 'fast_sqeuclidean') is unstable in this configuration"
        )

    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    # Haversine distance only accepts 2D data
    if metric == "haversine":
        X = X[:, :2]
        Y = Y[:, :2]

    argkmin_reduction = ArgKmin.get_for(
        X,
        Y,
        k=k,
        metric=metric,
        metric_kwargs=get_dummy_metric_kwargs(metric, n_features),
        # To be sure to use parallelization
        chunk_size=n_samples // 4,
    )

    dist_par_X, indices_par_X = argkmin_reduction.compute(
        strategy="parallel_on_X", return_distance=True
    )

    dist_par_Y, indices_par_Y = argkmin_reduction.compute(
        strategy="parallel_on_Y", return_distance=True
    )

    assert_argkmin_results_equality(
        dist_par_X, dist_par_Y, indices_par_X, indices_par_Y
    )


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("metric", RadiusNeighborhood.valid_metrics())
def test_radius_neighborhood_strategies_consistency(
    seed,
    n_samples,
    metric,
    n_features=10,
    radius=10.0,
    dtype=np.float64,
):
    # RadiusNeighborhood results obtained using both parallelization strategies
    # must be identical
    if _in_unstable_openblas_configuration() and metric == "fast_sqeuclidean":
        pytest.xfail(
            "OpenBLAS (used for 'fast_sqeuclidean') is unstable in this configuration"
        )

    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    # Haversine distance only accepts 2D data
    if metric == "haversine":
        X = X[:, :2]
        Y = Y[:, :2]

    radius_neigh_reduction = RadiusNeighborhood.get_for(
        X,
        Y,
        # Scaling the radius with the dimensions
        radius=radius ** np.log(n_features),
        metric=metric,
        metric_kwargs=get_dummy_metric_kwargs(metric, n_features),
        # To be sure to use parallelization
        chunk_size=n_samples // 4,
    )

    dist_par_X, indices_par_X = radius_neigh_reduction.compute(
        strategy="parallel_on_X", return_distance=True
    )

    dist_par_Y, indices_par_Y = radius_neigh_reduction.compute(
        strategy="parallel_on_Y", return_distance=True
    )

    assert_radius_neighborhood_results_equality(
        dist_par_X, dist_par_Y, indices_par_X, indices_par_Y
    )


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("n_features", [5, 10, 100])
@pytest.mark.parametrize("k, radius", [(50, 100)])
def test_fast_sqeuclidean_correctness(
    seed,
    n_samples,
    n_features,
    k,
    radius,
    dtype=np.float64,
):
    # The fast squared euclidean strategy must return results
    # that are close to the ones obtained with the euclidean distance
    if n_samples < k:
        pytest.skip(
            f"Skipping as n_samples (={n_samples}) < k (={k})",
            allow_module_level=True,
        )

    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    eucl_dist, eucl_indices = ArgKmin.get_for(X, Y, k, metric="euclidean").compute(
        return_distance=True
    )
    fse_dist, fse_indices = ArgKmin.get_for(X, Y, k, metric="fast_sqeuclidean").compute(
        return_distance=True
    )

    assert_argkmin_results_equality(eucl_dist, fse_dist, eucl_indices, fse_indices)

    eucl_dist, eucl_indices = RadiusNeighborhood.get_for(
        X, Y, radius, metric="euclidean"
    ).compute(return_distance=True)
    fse_dist, fse_indices = RadiusNeighborhood.get_for(
        X, Y, radius, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    assert_radius_neighborhood_results_equality(
        eucl_dist, fse_dist, eucl_indices, fse_indices
    )


@fails_if_unstable_openblas
@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3]])
@pytest.mark.parametrize("n_features", [5, 10, 100])
@pytest.mark.parametrize("k", [1, 10, 100])
@pytest.mark.parametrize("translation", [10 ** i for i in [4]])
def test_fast_sqeuclidean_translation_invariance(
    seed,
    n_samples,
    n_features,
    k,
    translation,
    dtype=np.float64,
):
    # The fast squared euclidean strategy should be translation invariant.
    if n_samples < k:
        pytest.skip(
            f"Skipping as n_samples (={n_samples}) < n_neighbors (={k})",
            allow_module_level=True,
        )

    rng = np.random.RandomState(seed)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    reference_dist, reference_indices = ArgKmin.get_for(
        X, Y, k, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    dist, indices = ArgKmin.get_for(
        X + translation, Y + translation, k, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    assert_argkmin_results_equality(reference_dist, dist, reference_indices, indices)
