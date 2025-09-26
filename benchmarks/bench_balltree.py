"""
================================================================
Benchmark: BallTree vs. scipy.spatial.cKDTree for k-NN Search
================================================================

This script compares the performance of sklearn.neighbors.BallTree with
scipy.spatial.cKDTree for a k-nearest neighbors (k-NN) query.

It measures and compares both the tree construction time and the query time
for each implementation.

Usage:
------
Run the script from the command line with optional parameters.

Example:
    # Run with 1000 samples and 100 features
    python bench_kdtree.py 1000 100

    # See all options
    python bench_kdtree.py --help
"""
import argparse
import time
import numpy as np
from sklearn.neighbors import BallTree
from scipy.spatial import cKDTree

def are_neighbors_equivalent(indices1: np.ndarray, indices2: np.ndarray) -> bool:
    """
    Checks if the two arrays of neighbor indices are equivalent.

    Since the order of neighbors at the same distance is not guaranteed,
    we sort each row before comparing.

    Parameters
    ----------
    indices1 : np.ndarray
        The neighbor indices from the first algorithm.
    indices2 : np.ndarray
        The neighbor indices from the second algorithm.

    Returns
    -------
    bool
        True if the sets of neighbors are the same for every point.
    """
    if indices1.shape != indices2.shape:
        return False
    # Sort the neighbors for each point to handle order differences
    indices1_sorted = np.sort(indices1, axis=1)
    indices2_sorted = np.sort(indices2, axis=1)
    return np.array_equal(indices1_sorted, indices2_sorted)


def run_benchmark(n_samples: int, n_features: int, leaf_size: int, k: int) -> None:
    """
    Generates data and runs the benchmark for BallTree and cKDTree.

    Parameters
    ----------
    n_samples : int
        Number of data points.
    n_features : int
        Number of dimensions for each point.
    leaf_size : int
        Leaf size parameter for both tree implementations.
    k : int
        Number of neighbors to query for.
    """
    # --- 1. Data Generation ---
    np.random.seed(42)  # for reproducibility
    X = np.random.rand(n_samples, n_features)

    print("-" * 60)
    print(f"Querying for {k} neighbors of {n_samples} points in {n_features} dimensions")
    print(f"(leaf_size = {leaf_size})")
    print("-" * 60)

    # --- 2. BallTree Benchmark ---
    t0 = time.time()
    bt = BallTree(X, leaf_size=leaf_size)
    construction_time_bt = time.time() - t0

    t0 = time.time()
    dist_bt, ind_bt = bt.query(X, k=k)
    query_time_bt = time.time() - t0
    total_time_bt = construction_time_bt + query_time_bt

    # --- 3. cKDTree Benchmark ---
    t0 = time.time()
    kdt = cKDTree(X, leafsize=leaf_size)
    construction_time_kdt = time.time() - t0

    t0 = time.time()
    # Note: cKDTree query k includes the point itself if it's in the dataset
    dist_kdt, ind_kdt = kdt.query(X, k=k)
    query_time_kdt = time.time() - t0
    total_time_kdt = construction_time_kdt + query_time_kdt

    # --- 4. Print Results ---
    print(f"{'Algorithm':<12} | {'Construction':<14} | {'Query':<14} | {'Total':<14}")
    print("-" * 60)
    print(f"{'BallTree':<12} | {construction_time_bt:<14.4f} | {query_time_bt:<14.4f} | {total_time_bt:<14.4f}")
    print(f"{'cKDTree':<12} | {construction_time_kdt:<14.4f} | {query_time_kdt:<14.4f} | {total_time_kdt:<14.4f}")
    print("-" * 60)

    # --- 5. Verify Results ---
    neighbors_match = are_neighbors_equivalent(ind_bt, ind_kdt)
    print(f"Neighbor sets match: {neighbors_match}")
    print("-" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark BallTree vs. cKDTree.")
    parser.add_argument(
        "n_samples", type=int, help="Number of samples to generate."
    )
    parser.add_argument(
        "n_features", type=int, help="Number of features (dimensions) for the data."
    )
    parser.add_argument(
        "--leaf_size", type=int, default=20, help="Leaf size for the trees."
    )
    parser.add_argument(
        "--k", type=int, default=20, dest="k_neighbors", help="Number of neighbors to query."
    )

    args = parser.parse_args()

    # Ensure k is not larger than the number of samples
    k = min(args.k_neighbors, args.n_samples)

    run_benchmark(
        n_samples=args.n_samples,
        n_features=args.n_features,
        leaf_size=args.leaf_size,
        k=k,
    )