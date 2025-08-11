"""
Test script to reproduce the HDBSCAN distance matrix modification issue #31907

This script demonstrates that HDBSCAN modifies the input precomputed distance matrix
when copy=False (the default behavior), which violates the principle that estimators
should not modify their input data unless explicitly requested.
"""

import numpy as np
from sklearn.cluster import HDBSCAN
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import make_blobs

# Generate sample data
X, y = make_blobs(n_samples=50, random_state=10)

# Create a precomputed distance matrix
D = euclidean_distances(X)
D_original = D.copy()

print("Testing HDBSCAN with precomputed distance matrix...")
print(f"Original matrix sum: {np.sum(D_original):.6f}")
print(f"Original matrix max: {np.max(D_original):.6f}")
print(f"Original matrix[0,1]: {D_original[0, 1]:.6f}")

# Test with copy=False (default behavior)
print("\n--- Testing with copy=False (default) ---")
D_test = D_original.copy()
print(f"Before HDBSCAN - matrix sum: {np.sum(D_test):.6f}")
print(f"Before HDBSCAN - matrix[0,1]: {D_test[0, 1]:.6f}")

clusterer = HDBSCAN(metric="precomputed", copy=False)
labels = clusterer.fit_predict(D_test)

print(f"After HDBSCAN - matrix sum: {np.sum(D_test):.6f}")
print(f"After HDBSCAN - matrix[0,1]: {D_test[0, 1]:.6f}")
print(f"Matrix modified: {not np.allclose(D_test, D_original)}")

# Test with copy=True
print("\n--- Testing with copy=True ---")
D_test = D_original.copy()
print(f"Before HDBSCAN - matrix sum: {np.sum(D_test):.6f}")
print(f"Before HDBSCAN - matrix[0,1]: {D_test[0, 1]:.6f}")

clusterer = HDBSCAN(metric="precomputed", copy=True)
labels = clusterer.fit_predict(D_test)

print(f"After HDBSCAN - matrix sum: {np.sum(D_test):.6f}")
print(f"After HDBSCAN - matrix[0,1]: {D_test[0, 1]:.6f}")
print(f"Matrix modified: {not np.allclose(D_test, D_original)}")

print("\n--- Analysis ---")
print("The issue is that when copy=False (default), HDBSCAN modifies the input matrix in-place")
print("This happens in two places:")
print("1. distance_matrix /= alpha  # Division by alpha parameter")
print("2. mutual_reachability_graph(distance_matrix, ...)  # In-place computation")
