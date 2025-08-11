"""
Test to reproduce HDBSCAN issue #31907 - modifying input precomputed distance matrix
"""
import numpy as np
import pytest
from sklearn.cluster import HDBSCAN
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import make_blobs
from sklearn.utils._testing import assert_allclose

def test_hdbscan_does_not_modify_precomputed_matrix_default():
    """Test that HDBSCAN with default copy=False does not modify input precomputed matrix.
    
    This test reproduces the issue #31907 where HDBSCAN modifies the input
    precomputed distance matrix when copy=False (default behavior).
    """
    # Generate sample data
    X, _ = make_blobs(n_samples=50, random_state=10)
    
    # Create precomputed distance matrix
    D = euclidean_distances(X)
    D_original = D.copy()
    
    # Test with default copy=False behavior
    # This should NOT modify the input matrix
    clusterer = HDBSCAN(metric="precomputed")  # copy=False by default
    _ = clusterer.fit_predict(D)
    
    # Assert that input matrix was not modified
    assert_allclose(D, D_original, err_msg="HDBSCAN modified input matrix with copy=False (default)")

def test_hdbscan_does_not_modify_precomputed_matrix_copy_false():
    """Test that HDBSCAN with explicit copy=False does not modify input precomputed matrix."""
    # Generate sample data
    X, _ = make_blobs(n_samples=50, random_state=10)
    
    # Create precomputed distance matrix
    D = euclidean_distances(X)
    D_original = D.copy()
    
    # Test with explicit copy=False
    # This should NOT modify the input matrix
    clusterer = HDBSCAN(metric="precomputed", copy=False)
    _ = clusterer.fit_predict(D)
    
    # Assert that input matrix was not modified
    assert_allclose(D, D_original, err_msg="HDBSCAN modified input matrix with copy=False")

def test_hdbscan_does_not_modify_precomputed_matrix_copy_true():
    """Test that HDBSCAN with copy=True does not modify input precomputed matrix."""
    # Generate sample data
    X, _ = make_blobs(n_samples=50, random_state=10)
    
    # Create precomputed distance matrix
    D = euclidean_distances(X)
    D_original = D.copy()
    
    # Test with copy=True
    # This should NOT modify the input matrix
    clusterer = HDBSCAN(metric="precomputed", copy=True)
    _ = clusterer.fit_predict(D)
    
    # Assert that input matrix was not modified
    assert_allclose(D, D_original, err_msg="HDBSCAN modified input matrix with copy=True")

if __name__ == "__main__":
    try:
        test_hdbscan_does_not_modify_precomputed_matrix_default()
        print("✓ Default behavior test passed")
    except AssertionError as e:
        print(f"✗ Default behavior test failed: {e}")
    
    try:
        test_hdbscan_does_not_modify_precomputed_matrix_copy_false()
        print("✓ copy=False test passed")
    except AssertionError as e:
        print(f"✗ copy=False test failed: {e}")
    
    try:
        test_hdbscan_does_not_modify_precomputed_matrix_copy_true()
        print("✓ copy=True test passed")
    except AssertionError as e:
        print(f"✗ copy=True test failed: {e}")
