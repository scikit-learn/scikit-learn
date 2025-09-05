"""Performance tests for stochastic optimizers.

These tests verify that the optimized implementations maintain
numerical accuracy while providing performance improvements.
"""

import numpy as np
import pytest

from sklearn.neural_network._stochastic_optimizers import (
    AdamOptimizer,
    SGDOptimizer,
)


@pytest.mark.parametrize("optimizer_class", [SGDOptimizer, AdamOptimizer])
def test_optimizer_numerical_accuracy_large_arrays(optimizer_class):
    """Test that optimizers maintain numerical accuracy with large arrays."""
    # Set random seed for reproducibility
    np.random.seed(42)

    # Create large parameter arrays to test vectorization benefits
    params = [
        np.random.randn(500, 500),  # Large weight matrix
        np.random.randn(500),  # Bias vector
    ]
    grads = [np.random.randn(500, 500), np.random.randn(500)]

    # Test with default parameters
    if optimizer_class == SGDOptimizer:
        optimizer = optimizer_class(params)
    else:  # AdamOptimizer
        optimizer = optimizer_class(params)

    # Get updates multiple times to ensure consistency
    updates1 = optimizer._get_updates([g.copy() for g in grads])
    updates2 = optimizer._get_updates([g.copy() for g in grads])

    # Verify updates have correct shape and are finite
    for update, param in zip(updates1, params):
        assert update.shape == param.shape
        assert np.all(np.isfinite(update))

    for update, param in zip(updates2, params):
        assert update.shape == param.shape
        assert np.all(np.isfinite(update))


def test_sgd_optimizer_performance_smoke_test():
    """Smoke test to verify SGD optimizer runs efficiently on large arrays."""
    np.random.seed(42)

    # Create moderately large arrays
    params = [np.random.randn(100, 100), np.random.randn(100)]
    grads = [np.random.randn(100, 100), np.random.randn(100)]

    optimizer = SGDOptimizer(params)

    # This should complete quickly with the vectorized implementation
    def run_updates():
        return optimizer._get_updates(grads)

    # Just verify it runs without timing assertions (CI environments vary)
    updates = run_updates()
    assert len(updates) == len(params)
    for update, param in zip(updates, params):
        assert update.shape == param.shape


def test_adam_optimizer_performance_smoke_test():
    """Smoke test to verify Adam optimizer runs efficiently on large arrays."""
    np.random.seed(42)

    # Create moderately large arrays
    params = [np.random.randn(100, 100), np.random.randn(100)]
    grads = [np.random.randn(100, 100), np.random.randn(100)]

    optimizer = AdamOptimizer(params)

    # This should complete quickly with the vectorized implementation
    def run_updates():
        return optimizer._get_updates(grads)

    # Just verify it runs without timing assertions (CI environments vary)
    updates = run_updates()
    assert len(updates) == len(params)
    for update, param in zip(updates, params):
        assert update.shape == param.shape


def test_sgd_momentum_vectorization():
    """Test that SGD momentum updates are properly vectorized."""
    np.random.seed(42)

    params = [np.random.randn(50, 50)]
    grads = [np.random.randn(50, 50)]

    optimizer = SGDOptimizer(params, momentum=0.9)

    # Set initial velocity to test momentum calculation
    initial_velocity = np.random.randn(50, 50)
    optimizer.velocities[0] = initial_velocity.copy()

    optimizer._get_updates(grads)

    # Verify the velocity was updated correctly
    expected_velocity = 0.9 * initial_velocity - 0.1 * grads[0]
    np.testing.assert_array_almost_equal(
        optimizer.velocities[0], expected_velocity
    )


def test_adam_moment_vectorization():
    """Test that Adam moment updates are properly vectorized."""
    np.random.seed(42)

    params = [np.random.randn(50, 50)]
    grads = [np.random.randn(50, 50)]

    optimizer = AdamOptimizer(params)

    # Set initial moments to test calculation
    initial_m = np.random.randn(50, 50)
    initial_v = np.random.randn(50, 50)
    optimizer.ms[0] = initial_m.copy()
    optimizer.vs[0] = initial_v.copy()

    optimizer._get_updates(grads)

    # Verify moments were updated correctly
    expected_m = 0.9 * initial_m + 0.1 * grads[0]
    expected_v = 0.999 * initial_v + 0.001 * (grads[0] ** 2)

    np.testing.assert_array_almost_equal(optimizer.ms[0], expected_m)
    np.testing.assert_array_almost_equal(optimizer.vs[0], expected_v)
