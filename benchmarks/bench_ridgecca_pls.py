"""
Benchmark of Ridge and SVD solvers for the CCA problem at different dimensions
"""
import gc
import sys
import numpy as np
from collections import defaultdict
from time import time
import matplotlib.pyplot as plt

from sklearn.datasets import make_regression
from sklearn.cross_decomposition import PLSSVD, RidgeCCA

plssvd = PLSSVD(n_components=1)
ridgecca = RidgeCCA(n_components=1, alpha_x=1.0, alpha_y=1.0)

def compute_bench(n_samples, features_range, n_repeats=5):
    results = defaultdict(lambda: [])

    for n_features in features_range:
        print("====================")
        print(f"Features: {n_features}")
        print("====================")

        for _ in range(n_repeats):
            dataset_kwargs = {
                "n_samples": n_samples,
                "n_features": n_features,
                "n_targets": n_features,
                "n_informative": n_features // 10,
                "bias": 0.0,
            }
            X, y = make_regression(**dataset_kwargs)

            for name, estimator in [("plssvd", plssvd), ("ridgecca", ridgecca)]:
                gc.collect()
                print(f"benchmarking {name}:", end="")
                sys.stdout.flush()
                tstart = time()
                estimator.fit(X, y)
                delta = time() - tstart
                print("%0.3fs" % delta)
                results[name].append(delta)

    # Convert lists to numpy arrays for easier manipulation
    for key in results:
        results[key] = np.array(results[key]).reshape(-1, n_repeats)

    return results

def plot_results(results, features_range):
    plt.figure(figsize=(10, 6))
    for method, times in results.items():
        mean_times = np.mean(times, axis=1)
        std_times = np.std(times, axis=1)
        plt.errorbar(features_range, mean_times, yerr=std_times, label=method, capsize=5)

    plt.xlabel('Number of Features')
    plt.ylabel('Time (seconds)')
    plt.title('Benchmarking Ridge and SVD Solvers for CCA')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    n_samples = 100
    features_range = [100, 200, 1000]
    results = compute_bench(n_samples, features_range)
    plot_results(results, features_range)
