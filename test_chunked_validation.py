#!/usr/bin/env python3
"""
Test script to demonstrate the new chunked validation functionality.

This script shows how the new chunk_size parameter in check_array and check_X_y
can help reduce memory usage when validating large datasets.
"""

import numpy as np
import sys
import os
import psutil
import time

# Add the sklearn directory to the path
sys.path.insert(0, '/Users/harshil/Documents/scikit-learn')

from sklearn.utils.validation import check_array, check_X_y
from sklearn.datasets import make_classification


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def test_chunked_validation():
    """Test the new chunked validation functionality."""
    print("Testing chunked validation functionality...")
    print("=" * 50)
    
    # Create a large dataset
    print("Creating large dataset (100,000 samples, 100 features)...")
    X, y = make_classification(
        n_samples=100000, 
        n_features=100, 
        n_informative=50,
        random_state=42
    )
    
    print(f"Dataset shape: {X.shape}")
    print(f"Initial memory usage: {get_memory_usage():.2f} MB")
    
    # Test 1: Standard validation (no chunking)
    print("\n1. Testing standard validation (no chunking)...")
    start_time = time.time()
    start_memory = get_memory_usage()
    
    try:
        X_checked_standard = check_array(
            X, 
            ensure_all_finite=True,
            ensure_non_negative=False
        )
        standard_time = time.time() - start_time
        standard_memory = get_memory_usage() - start_memory
        
        print(f"   ✓ Standard validation completed")
        print(f"   Time: {standard_time:.3f} seconds")
        print(f"   Memory increase: {standard_memory:.2f} MB")
        
    except Exception as e:
        print(f"   ✗ Standard validation failed: {e}")
        return False
    
    # Test 2: Chunked validation
    print("\n2. Testing chunked validation (chunk_size=10000)...")
    start_time = time.time()
    start_memory = get_memory_usage()
    
    try:
        X_checked_chunked = check_array(
            X, 
            ensure_all_finite=True,
            ensure_non_negative=False,
            chunk_size=10000
        )
        chunked_time = time.time() - start_time
        chunked_memory = get_memory_usage() - start_memory
        
        print(f"   ✓ Chunked validation completed")
        print(f"   Time: {chunked_time:.3f} seconds")
        print(f"   Memory increase: {chunked_memory:.2f} MB")
        
    except Exception as e:
        print(f"   ✗ Chunked validation failed: {e}")
        return False
    
    # Test 3: Verify results are identical
    print("\n3. Verifying results are identical...")
    if np.array_equal(X_checked_standard, X_checked_chunked):
        print("   ✓ Results are identical")
    else:
        print("   ✗ Results differ!")
        return False
    
    # Test 4: Test with check_X_y
    print("\n4. Testing chunked validation with check_X_y...")
    try:
        X_checked, y_checked = check_X_y(
            X, y,
            ensure_all_finite=True,
            chunk_size=10000
        )
        print(f"   ✓ check_X_y with chunking completed")
        print(f"   X shape: {X_checked.shape}, y shape: {y_checked.shape}")
    except Exception as e:
        print(f"   ✗ check_X_y with chunking failed: {e}")
        return False
    
    # Test 5: Test with non-negative validation
    print("\n5. Testing chunked non-negative validation...")
    # Create a dataset with some negative values to test the validation
    X_neg = X.copy()
    X_neg[0, 0] = -1.0  # Add a negative value
    
    try:
        # This should raise an error
        check_array(X_neg, ensure_non_negative=True, chunk_size=10000)
        print("   ✗ Should have raised an error for negative values")
        return False
    except ValueError as e:
        print(f"   ✓ Correctly caught negative value: {str(e)[:50]}...")
    
    print("\n" + "=" * 50)
    print("All tests passed! ✓")
    print("\nPerformance comparison:")
    print(f"  Standard validation: {standard_time:.3f}s, {standard_memory:.2f} MB")
    print(f"  Chunked validation:  {chunked_time:.3f}s, {chunked_memory:.2f} MB")
    
    if chunked_memory < standard_memory:
        print(f"  Memory reduction: {standard_memory - chunked_memory:.2f} MB")
    else:
        print("  Note: Memory usage may vary based on system state")
    
    return True


def test_edge_cases():
    """Test edge cases for chunked validation."""
    print("\nTesting edge cases...")
    print("=" * 30)
    
    # Test with small array (should not use chunking)
    print("1. Testing with small array...")
    X_small = np.random.rand(100, 10)
    try:
        X_checked = check_array(X_small, chunk_size=50)
        print("   ✓ Small array handled correctly")
    except Exception as e:
        print(f"   ✗ Small array failed: {e}")
        return False
    
    # Test with sparse matrix (should not use chunking)
    print("2. Testing with sparse matrix...")
    from scipy.sparse import csr_matrix
    X_sparse = csr_matrix(X_small)
    try:
        X_checked = check_array(X_sparse, chunk_size=50)
        print("   ✓ Sparse matrix handled correctly")
    except Exception as e:
        print(f"   ✗ Sparse matrix failed: {e}")
        return False
    
    # Test with invalid chunk_size
    print("3. Testing with invalid chunk_size...")
    try:
        X_checked = check_array(X_small, chunk_size=0)
        print("   ✓ Zero chunk_size handled correctly")
    except Exception as e:
        print(f"   ✗ Zero chunk_size failed: {e}")
        return False
    
    print("   ✓ All edge cases passed")
    return True


if __name__ == "__main__":
    print("Chunked Validation Test Suite")
    print("=" * 50)
    print("This script demonstrates the new chunked validation functionality")
    print("that helps reduce memory usage when validating large datasets.")
    print()
    
    success = True
    
    # Run main tests
    if not test_chunked_validation():
        success = False
    
    # Run edge case tests
    if not test_edge_cases():
        success = False
    
    if success:
        print("\n🎉 All tests passed! The chunked validation feature is working correctly.")
        print("\nThis enhancement will help users with large datasets avoid memory issues")
        print("during data validation by processing data in configurable chunks.")
    else:
        print("\n❌ Some tests failed. Please check the implementation.")
        sys.exit(1)
