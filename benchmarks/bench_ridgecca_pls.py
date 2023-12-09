import gc
import sys
from collections import defaultdict
from time import time

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cross_decomposition import PLSSVD, RidgeCCA, CCA
from sklearn.datasets import make_regression

# Initialize models
cca = CCA(n_components=1)
ridgecca_cca = RidgeCCA(n_components=1, alpha_x=0.0, alpha_y=0.0)
plssvd = PLSSVD(n_components=1)
ridgecca_pls = RidgeCCA(n_components=1, alpha_x=1.0, alpha_y=1.0)


def compute_bench(n_samples, features_range, methods, n_repeats=5):
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
                "n_informative": n_features,
            }
            X, y = make_regression(**dataset_kwargs)

            for name, estimator in methods:
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


def plot_results(results, features_range, title):
    plt.figure(figsize=(10, 6))
    markers = ['o', 's', 'D', '^',
               'x']  # You can add more markers if you have more methods
    for (method, times), marker in zip(results.items(), markers):
        mean_times = np.mean(times, axis=1)
        std_times = np.std(times, axis=1)
        plt.errorbar(
            features_range, mean_times, yerr=std_times, label=method, capsize=5,
            marker=marker
        )

    plt.xlabel("Number of Features", fontsize=14)
    plt.ylabel("Time (seconds)", fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    n_samples=500
    cca_features_range = [100, 200, 400]
    pls_features_range = [100, 200, 500, 1000, 2000]

    # Compute benchmarks for CCA
    cca_methods = [("cca", cca), ("ridgecca_cca", ridgecca_cca)]
    cca_results = compute_bench(n_samples, cca_features_range, cca_methods)
    plot_results(cca_results, cca_features_range, "Benchmarking CCA Methods")

    # Compute benchmarks for PLS
    pls_methods = [("plssvd", plssvd), ("ridgecca_pls", ridgecca_pls)]
    pls_results = compute_bench(n_samples, pls_features_range, pls_methods)
    plot_results(pls_results, pls_features_range, "Benchmarking PLS Methods")
