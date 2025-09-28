#!/usr/bin/env python3
"""
Example demonstrating the new chunked validation functionality.

This example shows how to use the new chunk_size parameter in check_array
and check_X_y to reduce memory usage when validating large datasets.
"""

import numpy as np
from sklearn.utils.validation import check_array, check_X_y
from sklearn.datasets import make_classification


def main():
    print("Chunked Validation Example")
    print("=" * 40)
    
    # Create a large dataset
    print("Creating large dataset (50,000 samples, 100 features)...")
    X, y = make_classification(
        n_samples=50000, 
        n_features=100, 
        n_informative=50,
        random_state=42
    )
    
    print(f"Dataset shape: {X.shape}")
    print(f"Memory usage: ~{X.nbytes / 1024 / 1024:.1f} MB")
    
    # Example 1: Standard validation
    print("\n1. Standard validation (no chunking):")
    try:
        X_checked = check_array(
            X, 
            ensure_all_finite=True,
            ensure_non_negative=True
        )
        print("   ✓ Validation completed successfully")
        print(f"   Result shape: {X_checked.shape}")
    except Exception as e:
        print(f"   ✗ Validation failed: {e}")
    
    # Example 2: Chunked validation
    print("\n2. Chunked validation (chunk_size=5000):")
    try:
        X_checked_chunked = check_array(
            X, 
            ensure_all_finite=True,
            ensure_non_negative=True,
            chunk_size=5000
        )
        print("   ✓ Chunked validation completed successfully")
        print(f"   Result shape: {X_checked_chunked.shape}")
        
        # Verify results are identical
        if np.array_equal(X_checked, X_checked_chunked):
            print("   ✓ Results are identical to standard validation")
        else:
            print("   ✗ Results differ from standard validation")
            
    except Exception as e:
        print(f"   ✗ Chunked validation failed: {e}")
    
    # Example 3: Using with check_X_y
    print("\n3. Chunked validation with check_X_y:")
    try:
        X_checked, y_checked = check_X_y(
            X, y,
            ensure_all_finite=True,
            chunk_size=5000
        )
        print("   ✓ check_X_y with chunking completed successfully")
        print(f"   X shape: {X_checked.shape}, y shape: {y_checked.shape}")
    except Exception as e:
        print(f"   ✗ check_X_y with chunking failed: {e}")
    
    # Example 4: Error handling with chunked validation
    print("\n4. Error handling with chunked validation:")
    X_with_neg = X.copy()
    X_with_neg[0, 0] = -1.0  # Add a negative value
    
    try:
        check_array(X_with_neg, ensure_non_negative=True, chunk_size=5000)
        print("   ✗ Should have raised an error for negative values")
    except ValueError as e:
        print(f"   ✓ Correctly caught negative value: {str(e)[:60]}...")
    
    print("\n" + "=" * 40)
    print("Example completed!")
    print("\nKey benefits of chunked validation:")
    print("- Reduces memory usage for large datasets")
    print("- Maintains identical results to standard validation")
    print("- Helps prevent memory errors with very large arrays")
    print("- Works with both check_array and check_X_y")


if __name__ == "__main__":
    main()
