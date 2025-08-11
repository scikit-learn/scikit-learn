"""
Simple demonstration of the HDBSCAN fix for issue #31907
This shows the difference between the old buggy behavior and the new fixed behavior.
"""

import numpy as np

def euclidean_distances(X):
    """Simple euclidean distance matrix computation"""
    n = X.shape[0]
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            distances[i, j] = np.sqrt(np.sum((X[i] - X[j])**2))
    return distances

def old_hdbscan_brute(X, alpha=1.0, copy=False):
    """Old implementation that has the bug"""
    # This is the problematic line from the original code
    distance_matrix = X.copy() if copy else X  # BUG: when copy=False, modifies X directly
    distance_matrix /= alpha  # This modifies the matrix in-place
    return distance_matrix

def fixed_hdbscan_brute(X, alpha=1.0, copy=False):
    """Fixed implementation"""
    # For precomputed distance matrices, always copy to avoid modifying input data.
    # This ensures scikit-learn's principle of not modifying input data is preserved.
    # The copy parameter is kept for backwards compatibility.
    distance_matrix = X.copy()  # FIXED: Always copy for precomputed matrices
    distance_matrix /= alpha  # This modifies the matrix in-place, but it's a copy
    return distance_matrix

def test_fix():
    """Test the fix"""
    print("Testing HDBSCAN fix for issue #31907\n")
    
    # Create sample data
    np.random.seed(10)
    X_features = np.random.randn(5, 2)
    D = euclidean_distances(X_features)
    D_original = D.copy()
    
    print("Original distance matrix:")
    print(f"Sum: {np.sum(D_original):.6f}")
    print(f"Element [0,1]: {D_original[0, 1]:.6f}\n")
    
    # Test old implementation (copy=False, default)
    print("=== OLD Implementation (copy=False) ===")
    D_test_old = D_original.copy()
    print("Before:", f"Sum: {np.sum(D_test_old):.6f}, Element [0,1]: {D_test_old[0, 1]:.6f}")
    
    old_hdbscan_brute(D_test_old, alpha=2.0, copy=False)
    print("After: ", f"Sum: {np.sum(D_test_old):.6f}, Element [0,1]: {D_test_old[0, 1]:.6f}")
    print(f"Matrix modified: {not np.allclose(D_test_old, D_original)}")
    
    # Test fixed implementation (copy=False, default)
    print("\n=== FIXED Implementation (copy=False) ===")
    D_test_fixed = D_original.copy()
    print("Before:", f"Sum: {np.sum(D_test_fixed):.6f}, Element [0,1]: {D_test_fixed[0, 1]:.6f}")
    
    fixed_hdbscan_brute(D_test_fixed, alpha=2.0, copy=False)
    print("After: ", f"Sum: {np.sum(D_test_fixed):.6f}, Element [0,1]: {D_test_fixed[0, 1]:.6f}")
    print(f"Matrix modified: {not np.allclose(D_test_fixed, D_original)}")
    
    # Results
    old_modified = not np.allclose(D_test_old, D_original)
    fixed_modified = not np.allclose(D_test_fixed, D_original)
    
    print("\n=== RESULTS ===")
    if old_modified and not fixed_modified:
        print("✅ SUCCESS: Fix works correctly!")
        print("   - Old code incorrectly modifies input (BUG)")
        print("   - Fixed code preserves input (CORRECT)")
    elif not old_modified:
        print("❌ Test setup issue: Old code should have modified the input")
    elif fixed_modified:
        print("❌ Fix didn't work: Fixed code still modifies input")
    else:
        print("❓ Unexpected result")
    
    print(f"\nOld implementation modified input: {old_modified}")
    print(f"Fixed implementation modified input: {fixed_modified}")
    
    return old_modified and not fixed_modified

if __name__ == "__main__":
    test_fix()
