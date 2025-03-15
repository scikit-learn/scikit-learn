"""Tests for TSNEParticleSwarmOptimizer."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

# No need for path manipulation since we're using absolute imports
from sklearn.datasets import load_iris
from sklearn.manifold._tsne_pso import TSNEPSO
from sklearn.preprocessing import StandardScaler


# Define a custom assert_no_warnings function since it's not available in this version
def assert_no_warnings(func, *args, **kw):
    """Assert that no warnings are raised during function execution."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        result = func(*args, **kw)
    assert len(record) == 0, f"Got warnings: {[str(w.message) for w in record]}"
    return result


def test_tsnepso_init_params():
    """Test TSNEPSO initialization parameters."""
    # Test default parameters
    tsne_pso = TSNEPSO()
    assert tsne_pso.n_components == 2
    assert tsne_pso.perplexity == 30.0
    assert tsne_pso.n_particles == 10
    assert tsne_pso.use_hybrid is True

    # Test custom parameters
    tsne_pso = TSNEPSO(
        n_components=3,
        perplexity=10.0,
        n_iter=500,
        n_particles=5,
        inertia_weight=0.3,
        cognitive_weight=0.8,
        social_weight=1.2,
        use_hybrid=False,
        degrees_of_freedom=0.5,
        init="random",
        verbose=1,
        random_state=42,
    )
    assert tsne_pso.n_components == 3
    assert tsne_pso.perplexity == 10.0
    assert tsne_pso.n_iter == 500
    assert tsne_pso.n_particles == 5
    assert tsne_pso.inertia_weight == 0.3
    assert tsne_pso.cognitive_weight == 0.8
    assert tsne_pso.social_weight == 1.2
    assert tsne_pso.use_hybrid is False
    assert tsne_pso.degrees_of_freedom == 0.5
    assert tsne_pso.init == "random"
    assert tsne_pso.verbose == 1
    assert tsne_pso.random_state == 42


def test_tsnepso_validation():
    """Test parameter validation in TSNEPSO."""
    # Test invalid n_components
    tsne_pso = TSNEPSO(n_components=0)
    with pytest.raises(ValueError, match="n_components.*range.*Got 0"):
        tsne_pso._validate_parameters()

    # Test invalid perplexity
    tsne_pso = TSNEPSO(perplexity=0)
    with pytest.raises(ValueError, match="perplexity must be greater than 0"):
        tsne_pso._validate_parameters()

    # Test invalid n_iter
    tsne_pso = TSNEPSO(n_iter=0)
    with pytest.raises(ValueError, match="n_iter.*range"):
        tsne_pso._validate_parameters()

    # Test invalid method
    tsne_pso = TSNEPSO(method="invalid")
    with pytest.raises(ValueError, match="must be .+ 'pso'"):
        tsne_pso._validate_parameters()

    # Test invalid early exaggeration
    tsne_pso = TSNEPSO(early_exaggeration=0)
    with pytest.raises(ValueError, match="early_exaggeration.*positive"):
        tsne_pso._validate_parameters()

    # Test invalid n_particles
    tsne_pso = TSNEPSO(n_particles=0)
    with pytest.raises(ValueError, match="n_particles.*range"):
        tsne_pso._validate_parameters()

    # Test invalid inertia_weight
    tsne_pso = TSNEPSO(inertia_weight=1.5)
    with pytest.raises(ValueError, match="inertia_weight.*range"):
        tsne_pso._validate_parameters()


def test_tsnepso_iris():
    """Test TSNEPSO on the Iris dataset."""
    # Load and scale the data
    iris = load_iris()
    X = StandardScaler().fit_transform(iris.data)

    # Fit TSNEPSO with reduced iterations for testing
    tsne_pso = TSNEPSO(
        n_components=2, perplexity=10.0, n_iter=50, n_particles=3, random_state=42
    )
    embedding = tsne_pso.fit_transform(X)

    # Check output shape
    assert embedding.shape == (X.shape[0], 2)

    # Check attributes
    assert hasattr(tsne_pso, "embedding_")
    assert hasattr(tsne_pso, "kl_divergence_")
    assert hasattr(tsne_pso, "n_iter_")
    assert hasattr(tsne_pso, "n_features_in_")

    # Check n_features_in_
    assert tsne_pso.n_features_in_ == X.shape[1]


def test_tsnepso_iris_init_array():
    """Test TSNEPSO with array initialization."""
    # Load and scale the data
    iris = load_iris()
    X = StandardScaler().fit_transform(iris.data)

    # Create an initial embedding
    init_embedding = np.random.RandomState(42).normal(0, 0.0001, (X.shape[0], 2))

    # Fit TSNEPSO with the initial embedding
    tsne_pso = TSNEPSO(
        n_components=2,
        perplexity=10.0,
        n_iter=50,
        n_particles=3,
        init=init_embedding,
        random_state=42,
    )
    embedding = tsne_pso.fit_transform(X)

    # Check output shape
    assert embedding.shape == (X.shape[0], 2)

    # Check attributes
    assert hasattr(tsne_pso, "embedding_")
    assert hasattr(tsne_pso, "kl_divergence_")

    # Check for wrong init shape
    wrong_shape_init = np.random.RandomState(42).normal(0, 0.0001, (X.shape[0], 3))
    tsne_pso = TSNEPSO(init=wrong_shape_init, n_components=2)
    with pytest.raises(ValueError, match="init.shape="):
        tsne_pso.fit_transform(X)


def test_tsnepso_transform_raises():
    """Test that transform raises NotImplementedError."""
    tsne_pso = TSNEPSO()
    with pytest.raises(NotImplementedError):
        tsne_pso.transform(np.random.random((10, 5)))


def test_tsnepso_precomputed():
    """Test TSNEPSO with precomputed distances."""
    # Load and scale the data
    iris = load_iris()
    X = StandardScaler().fit_transform(iris.data)

    # Compute pairwise distances
    from sklearn.metrics import pairwise_distances

    distances = pairwise_distances(X, metric="euclidean", squared=True)

    # Fit TSNEPSO with precomputed distances
    tsne_pso = TSNEPSO(
        n_components=2,
        perplexity=10.0,
        n_iter=50,
        n_particles=3,
        metric="precomputed",
        random_state=42,
    )
    embedding = tsne_pso.fit_transform(distances)

    # Check output shape
    assert embedding.shape == (X.shape[0], 2)

    # Check for non-square precomputed distances
    non_square = np.random.random((10, 8))
    tsne_pso = TSNEPSO(metric="precomputed")
    with pytest.raises(ValueError, match="X should be a square distance matrix"):
        tsne_pso.fit_transform(non_square)

    # Check for negative values in precomputed distances
    negative_dists = np.random.random((10, 10))
    negative_dists[0, 1] = -1
    tsne_pso = TSNEPSO(metric="precomputed")
    with pytest.raises(
        ValueError, match="Precomputed distance contains negative values"
    ):
        tsne_pso.fit_transform(negative_dists)


def test_tsnepso_perplexity_warning():
    """Test warning when perplexity is too high for the number of samples."""
    # Create a small dataset
    X = np.random.random((10, 5))

    # Set perplexity too high
    tsne_pso = TSNEPSO(perplexity=5)

    # Check that warning is raised
    with pytest.warns(
        UserWarning, match="Perplexity is too large for the number of samples"
    ):
        tsne_pso.fit_transform(X)


def test_tsnepso_random_state():
    """Test that random_state ensures reproducible results."""
    # Load and scale the data
    iris = load_iris()
    X = StandardScaler().fit_transform(iris.data)

    # Fit with fixed random state
    tsne_pso1 = TSNEPSO(
        n_components=2, perplexity=10.0, n_iter=50, n_particles=3, random_state=42
    )
    embedding1 = tsne_pso1.fit_transform(X)

    # Fit again with the same random state
    tsne_pso2 = TSNEPSO(
        n_components=2, perplexity=10.0, n_iter=50, n_particles=3, random_state=42
    )
    embedding2 = tsne_pso2.fit_transform(X)

    # Results should be exactly the same
    assert_allclose(embedding1, embedding2)
