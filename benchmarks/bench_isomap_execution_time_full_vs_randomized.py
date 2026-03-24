"""
======================================================================
Isomap Solvers Benchmark: Execution Time vs Number of Samples
======================================================================

This benchmark demonstrates how the choice of eigen_solver in Isomap
can significantly affect computation time, especially as the dataset
size increases.

Description:
------------
Synthetic datasets are generated using `make_classification` with a
fixed number of features. The number of samples is
varied from 1000 to 4000.

For each setting, Isomap is applied using two different solvers:
- 'auto' (full eigendecomposition)
- 'randomized_value'

The execution time of each solver is measured for each number of
samples, and the average time over multiple runs (default: 3) is
plotted.

What you can observe:
---------------------
If n_components < 10, the randomized and auto solvers produce similar
results (in this case, the arpack solver is selected).
However, when n_components > 10, the randomized solver becomes significantly
faster, especially as the number of samples increases.

"""

import time

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification
from sklearn.manifold import Isomap

# 1 - Experiment Setup
# -- -- -- -- -- -- -- -- -- -
n_samples_list = [1000, 2000, 3000, 4000]
n_neighbors = 30
n_components_list = [2, 10]
n_features = 100
n_iter = 3  # Number of repetitions for averaging execution time

# Store timings for each value of n_components
timing_all = {}

for n_components in n_components_list:
    # Create containers for timing results
    timing = {
        "auto": np.zeros((len(n_samples_list), n_iter)),
        "randomized_value": np.zeros((len(n_samples_list), n_iter)),
    }

    for j, n in enumerate(n_samples_list):
        # Generate synthetic classification dataset
        X, _ = make_classification(
            n_samples=n,
            n_features=n_features,
            n_redundant=0,
            n_clusters_per_class=1,
            n_classes=1,
            random_state=42,
        )

        # Evaluate both solvers for multiple repetitions
        for solver in ["auto", "randomized_value"]:
            for i in range(n_iter):
                model = Isomap(
                    n_neighbors=n_neighbors,
                    n_components=n_components,
                    eigen_solver=solver,
                )
                start = time.perf_counter()
                model.fit(X)
                elapsed = time.perf_counter() - start
                timing[solver][j, i] = elapsed

    timing_all[n_components] = timing

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, n_components in enumerate(n_components_list):
    ax = axes[idx]
    timing = timing_all[n_components]
    avg_full = timing["auto"].mean(axis=1)
    std_full = timing["auto"].std(axis=1)
    avg_rand = timing["randomized_value"].mean(axis=1)
    std_rand = timing["randomized_value"].std(axis=1)

    ax.errorbar(
        n_samples_list,
        avg_full,
        yerr=std_full,
        label="Isomap (full)",
        marker="o",
        linestyle="-",
    )
    ax.errorbar(
        n_samples_list,
        avg_rand,
        yerr=std_rand,
        label="Isomap (randomized)",
        marker="x",
        linestyle="--",
    )
    ax.set_xlabel("Number of Samples")
    ax.set_ylabel("Execution Time (seconds)")
    ax.set_title(f"Isomap Execution Time (n_components = {n_components})")
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()
