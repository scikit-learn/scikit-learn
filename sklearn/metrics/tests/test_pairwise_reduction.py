import itertools

import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_allclose
from scipy.sparse import (
    lil_matrix,
    csc_matrix,
    dok_matrix,
    dia_matrix,
    coo_matrix,
    bsr_matrix,
)

from sklearn.metrics._parallel_reductions import (
    ArgKmin,
    RadiusNeighborhood,
    FastSquaredEuclideanArgKmin,
    SparseSparseDatasetsPair,
    DenseDenseDatasetsPair,
    DenseSparseDatasetsPair,
    SparseDenseDatasetsPair,
    FastSquaredEuclideanRadiusNeighborhood,
)


def assert_radius_neighborhood_equality(ref_dist, dist, ref_indices, indices):
    # We get arrays of arrays and we need to check for individual pairs
    for i in range(ref_dist.shape[0]):
        assert_array_equal(
            ref_dist[i],
            dist[i],
            err_msg=f"Query vector #{i} has different neighbors' distances",
        )
        assert_array_equal(
            ref_indices[i],
            indices[i],
            err_msg=f"Query vector #{i} has different neighbors' indices",
        )


def test_argkmin_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    k = 5
    metric = "euclidean"

    with pytest.raises(ValueError, match="`k`= -1, must be >= 1."):
        ArgKmin.get_for(X=X, Y=Y, k=-1, metric=metric)

    with pytest.raises(ValueError, match="`k`= 0, must be >= 1."):
        ArgKmin.get_for(X=X, Y=Y, k=0.1, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        ArgKmin.get_for(X=X, Y=Y, k=k, metric="wrong metric")

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        ArgKmin.get_for(X=[1, 2], Y=Y, k=k, metric=metric)

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        ArgKmin.get_for(X=X, Y=[1, 2], k=k, metric=metric)


def test_radius_neighborhood_factory_method_wrong_usages():
    rng = np.random.RandomState(1)
    X = rng.rand(100, 10)
    Y = rng.rand(100, 10)
    radius = 5
    metric = "euclidean"

    with pytest.raises(ValueError, match="`radius`= -1.0, must be >= 0."):
        RadiusNeighborhood.get_for(X=X, Y=Y, radius=-1, metric=metric)

    with pytest.raises(ValueError, match="Unrecognized metric"):
        RadiusNeighborhood.get_for(X=X, Y=Y, radius=radius, metric="wrong metric")

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        RadiusNeighborhood.get_for(X=[1, 2], Y=Y, radius=radius, metric=metric)

    with pytest.raises(ValueError, match="Expected 2D array, got 1D array instead"):
        RadiusNeighborhood.get_for(X=X, Y=[1, 2], radius=radius, metric=metric)


@pytest.mark.parametrize(
    "PairwiseDistancesReduction, FastSquaredPairwiseDistancesReduction",
    [
        (ArgKmin, FastSquaredEuclideanArgKmin),
        (RadiusNeighborhood, FastSquaredEuclideanRadiusNeighborhood),
    ],
)
def test_paiwise_disances_reduction_factory_method(
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


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("k", [1, 10, 100])
@pytest.mark.parametrize("chunk_size", [512, 1024, 1337, 19301])
@pytest.mark.parametrize("metric", ["euclidean", "fast_sqeuclidean"])
def test_argkmin_chunk_size_agnosticism(
    n_samples, k, chunk_size, metric, n_features=100, dtype=np.float64
):
    # ArgKmin results should not depend on the chunk size
    rng = np.random.RandomState(1)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    ref_dist, ref_indices = ArgKmin.get_for(X, Y, k=k, metric="euclidean").compute(
        return_distance=True
    )

    dist, indices = ArgKmin.get_for(
        X, Y, k=k, metric=metric, chunk_size=chunk_size
    ).compute(return_distance=True)

    assert_array_equal(ref_dist, dist)
    assert_array_equal(ref_indices, indices)


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("radius", [1, 10, 100])
@pytest.mark.parametrize("chunk_size", [512, 1024, 1337, 19301])
@pytest.mark.parametrize("metric", ["euclidean", "fast_sqeuclidean"])
def test_radius_neighborhood_chunk_size_agnosticism(
    n_samples, radius, chunk_size, metric, n_features=100, dtype=np.float64
):
    # RadiusNeighborhood results should not depend on the chunk size
    rng = np.random.RandomState(1)
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

    assert_radius_neighborhood_equality(ref_dist, dist, ref_indices, indices)


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("n_features", [5, 100, 500])
@pytest.mark.parametrize("k", [1, 10, 100])
@pytest.mark.parametrize("metric", ["euclidean", "fast_sqeuclidean"])
def test_argkmin_strategies_consistency(
    n_samples,
    n_features,
    k,
    metric,
    dtype=np.float64,
):
    # ArgKmin results obtained using both parallelization strategies
    # must be identical

    rng = np.random.RandomState(1)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    argkmin_reduction = ArgKmin.get_for(X, Y, k=k, metric=metric)

    dist_par_X, indices_par_X = argkmin_reduction.compute(
        strategy="parallel_on_X", return_distance=True
    )

    dist_par_Y, indices_par_Y = argkmin_reduction.compute(
        strategy="parallel_on_Y", return_distance=True
    )

    assert_array_equal(dist_par_X, dist_par_Y)
    assert_array_equal(indices_par_X, indices_par_Y)


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("n_features", [5, 100, 500])
@pytest.mark.parametrize("radius", [1, 10, 100])
@pytest.mark.parametrize("metric", ["euclidean", "fast_sqeuclidean"])
def test_radius_neighborhood_strategies_consistency(
    n_samples,
    n_features,
    radius,
    metric,
    dtype=np.float64,
):
    # RadiusNeighborhood results obtained using both parallelization strategies
    # must be identical

    rng = np.random.RandomState(1)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    radius_neigh_reduction = RadiusNeighborhood.get_for(
        X,
        Y,
        # Scaling the radius with the dimensions
        radius=radius ** np.log(n_features),
        metric=metric,
    )

    dist_par_X, indices_par_X = radius_neigh_reduction.compute(
        strategy="parallel_on_X", return_distance=True
    )

    dist_par_Y, indices_par_Y = radius_neigh_reduction.compute(
        strategy="parallel_on_Y", return_distance=True
    )

    assert_radius_neighborhood_equality(
        dist_par_X, dist_par_Y, indices_par_X, indices_par_Y
    )


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("n_features", [5, 10, 100])
@pytest.mark.parametrize("sample_imbalance", [10, 2, 1, 0.5])
@pytest.mark.parametrize("k, radius", [(1, 0), (10, 1), (100, 10), (1000, 100)])
def test_fast_sqeuclidean_correctness(
    n_samples,
    n_features,
    sample_imbalance,
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

    rng = np.random.RandomState(1)

    spread = 100
    X = (
        rng.rand(int(n_samples * n_features / sample_imbalance))
        .astype(dtype)
        .reshape((-1, n_features))
        * spread
    )
    Y = (
        rng.rand(int(n_samples * n_features)).astype(dtype).reshape((-1, n_features))
        * spread
    )

    eucl_dist, eucl_indices = ArgKmin.get_for(X, Y, k, metric="euclidean").compute(
        return_distance=True
    )
    fse_dist, fse_indices = ArgKmin.get_for(X, Y, k, metric="fast_sqeuclidean").compute(
        return_distance=True
    )

    assert_array_equal(eucl_dist, fse_dist)
    assert_array_equal(eucl_indices, fse_indices)

    eucl_dist, eucl_indices = RadiusNeighborhood.get_for(
        X, Y, radius, metric="euclidean"
    ).compute(return_distance=True)
    fse_dist, fse_indices = RadiusNeighborhood.get_for(
        X, Y, radius, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    assert_radius_neighborhood_equality(eucl_dist, fse_dist, eucl_indices, fse_indices)


@pytest.mark.parametrize("n_samples", [10 ** i for i in [2, 3, 4]])
@pytest.mark.parametrize("n_features", [5, 10, 100, 500])
@pytest.mark.parametrize("k", [1, 10, 100, 1000])
@pytest.mark.parametrize("translation", [10 ** i for i in [2, 3, 4, 5, 6, 7]])
@pytest.mark.skip(
    reason=(
        "Long test, translation invariance should have its own study: skipping for now"
    )
)
def test_fast_sqeuclidean_translation_invariance(
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

    rng = np.random.RandomState(1)
    spread = 100
    X = rng.rand(n_samples, n_features).astype(dtype) * spread
    Y = rng.rand(n_samples, n_features).astype(dtype) * spread

    reference_dist, reference_indices = ArgKmin.get_for(
        X, Y, k, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    dist, indices = ArgKmin.get_for(
        X + translation, X + translation, k, metric="fast_sqeuclidean"
    ).compute(return_distance=True)

    assert_array_equal(reference_indices, indices)
    assert_allclose(reference_dist, dist)
