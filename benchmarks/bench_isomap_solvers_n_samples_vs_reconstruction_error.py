"""
========================================================================
Benchmark: Isomap Reconstruction Error - Standard vs. Randomized Solver
========================================================================

This benchmark illustrates how the number of samples impacts the quality 
of the Isomap embedding, using reconstruction error as a metric.

Description:
------------
We generate synthetic 2D non-linear data (two concentric circles) with 
varying numbers of samples. For each subset, we compare the reconstruction 
error of two Isomap solvers:

- The `auto` solver (standard dense or arpack, selected automatically).
- The  `randomized_value` solver .

What you can observe:
---------------------
- The difference in performance between the two solvers.

Further exploration:
---------------------
- Modify the number of neighbors or iterations.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import Isomap
from sklearn.datasets import make_circles

# 1- Experiment Configuration
# ---------------------------
min_n_samples, max_n_samples = 100, 4000
n_samples_grid_size = 4  # Number of sample sizes to test

n_samples_range = [
    int(min_n_samples + np.floor((x / (n_samples_grid_size - 1)) * (max_n_samples - min_n_samples)))
    for x in range(0, n_samples_grid_size)
]

n_components = 2
n_iter = 3  # Number of repetitions per sample size
include_arpack = False  # Reserved for further testing

# 2- Data Generation
# ------------------
n_features = 2
X_full, y_full = make_circles(n_samples=max_n_samples, factor=0.3, noise=0.05, random_state=0)

# 3- Benchmark Execution
# ----------------------
errors_randomized = []
errors_full = []

for n_samples in n_samples_range:
    X, y = X_full[:n_samples], y_full[:n_samples]
    print(f"Computing for n_samples = {n_samples}")

    # Instantiate Isomap solvers
    isomap_randomized = Isomap(n_neighbors=50, n_components=n_components, eigen_solver='randomized_value')
    isomap_auto = Isomap(n_neighbors=50, n_components=n_components, eigen_solver='auto')

    # Fit and record reconstruction error
    isomap_randomized.fit(X)
    err_rand = isomap_randomized.reconstruction_error()
    errors_randomized.append(err_rand)

    isomap_auto.fit(X)
    err_auto = isomap_auto.reconstruction_error()
    errors_full.append(err_auto)

# 4- Results Visualization
# ------------------------
plt.figure(figsize=(10, 6))
plt.scatter(n_samples_range, errors_full, label='Isomap (auto)', color='b', marker='*')
plt.scatter(n_samples_range, errors_randomized, label='Isomap (randomized)', color='r', marker='x')

plt.title('Isomap Reconstruction Error vs. Number of Samples')
plt.xlabel('Number of Samples')
plt.ylabel('Reconstruction Error')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
