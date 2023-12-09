"""
Benchmarking CCA/PLS Methods with Constant Sample Size

In this benchmark we show that PCA-CCA and PCA-PLS (i.e. solving the CCA and PLS
problems in the PCA space) are faster than the standard CCA and PLS algorithms,
respectively for high-dimensional data.
The RidgeCCA estimator uses PCA behind the scenes and thus can be applied by setting
the alpha_x and alpha_y parameters to 0.0 (CCA) or 1.0 (PLS).

PCA-CCA consistently outperforms the standard CCA algorithm in terms of computation
time across the entire range of feature sizes. Notably, the standard CCA method is
constrained to datasets where the number of features is less than the number of
samples (p < n), which limits the maximum feature size to 400 in our tests. This
limitation is due to the requirement of CCA for the inversion of the covariance
matrix, which must be full rank and thus invertible.

As observed in the plot, the PLSSVD method initially performs faster than the
pca_pls method. However, as the number of features increases, we notice a
'crossover' point: the computation time for PLSSVD grows more rapidly than for
pca_pls. Beyond this crossover, pca_pls becomes the more efficient
algorithm in terms of computation time. This suggests that pca_pls may be better
suited for high-dimensional data scenarios typically encountered in modern datasets.

It's important to note that for CCA problems, which are closely related to PLS,
 the Ridge CCA method demonstrates faster computation times across all feature sizes
tested. However, unlike PLS, CCA is only uniquely defined when the number of features
(p) is less than the number of samples (n), which constrained our tests to a maximum
of 400 features for the CCA benchmark.
"""

import gc
import sys
from collections import defaultdict
from time import time

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cross_decomposition import CCA, PLSSVD, RidgeCCA
from sklearn.datasets import make_regression

# Initialize models
cca = CCA(n_components=1)
pca_cca = RidgeCCA(n_components=1, alpha_x=0.0, alpha_y=0.0)
plssvd = PLSSVD(n_components=1)
pca_pls = RidgeCCA(n_components=1, alpha_x=1.0, alpha_y=1.0)


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
    markers = [
        "o",
        "s",
        "D",
        "^",
        "x",
    ]  # You can add more markers if you have more methods
    for (method, times), marker in zip(results.items(), markers):
        mean_times = np.mean(times, axis=1)
        std_times = np.std(times, axis=1)
        plt.errorbar(
            features_range,
            mean_times,
            yerr=std_times,
            label=method,
            capsize=5,
            marker=marker,
        )

    plt.xlabel("Number of Features", fontsize=14)
    plt.ylabel("Time (seconds)", fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend()
    plt.grid(True)


if __name__ == "__main__":
    n_samples = 500
    cca_features_range = [100, 200, 400]
    pls_features_range = [100, 200, 500, 1000, 2000]

    # Compute benchmarks for CCA
    cca_methods = [("cca", cca), ("pca-cca", pca_cca)]
    cca_results = compute_bench(n_samples, cca_features_range, cca_methods)
    plot_results(cca_results, cca_features_range, "Benchmarking CCA Methods")
    plt.show()

    # Compute benchmarks for PLS
    pls_methods = [("plssvd", plssvd), ("pca-pls", pca_pls)]
    pls_results = compute_bench(n_samples, pls_features_range, pls_methods)
    plot_results(pls_results, pls_features_range, "Benchmarking PLS Methods")
    plt.show()
