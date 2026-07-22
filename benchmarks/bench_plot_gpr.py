"""
Timing Benchmark Test for GaussianProcessRegressor.

This Benchmark Test is used to evaluate the effect of different values
of n_jobs parameter on the performance of the regression algorithm.
"""

import time

import matplotlib.pyplot as plt
import pandas as pd
from joblib import parallel_backend

from sklearn.datasets import make_regression
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C

fixed_kernel = RBF(length_scale=1.0, length_scale_bounds="fixed")
kernel = C(0.1, (1e-2, 1e2)) * RBF(
    length_scale=1.0, length_scale_bounds=(1e-3, 1e3)
) + C(1e-5, (1e-5, 1e2))

max_n_restarts = 5
n_trials = 5

X, y = make_regression(n_samples=500, n_features=5, n_informative=5, noise=1.0)
run_data = []

for n_restarts in range(max_n_restarts + 1):
    for trial in range(n_trials):
        start = time.perf_counter()
        with parallel_backend("sequential", n_jobs=1):
            GaussianProcessRegressor(
                kernel=kernel, n_restarts_optimizer=n_restarts, random_state=42
            ).fit(X, y)
        run_data.append(
            {
                "backend": "sequential",
                "n_restarts": n_restarts,
                "execution_time": time.perf_counter() - start,
            }
        )

    for trial in range(n_trials):
        start = time.perf_counter()
        with parallel_backend("loky", n_jobs=-1):
            GaussianProcessRegressor(
                kernel=kernel, n_restarts_optimizer=n_restarts, random_state=42
            ).fit(X, y)
        run_data.append(
            {
                "backend": "parallel (loky)",
                "n_restarts": n_restarts,
                "execution_time": time.perf_counter() - start,
            }
        )

df = pd.DataFrame(run_data)

grouped = df.groupby(["backend", "n_restarts"])["execution_time"]
stats = grouped.agg(["mean", "std", "count"]).reset_index()

fig, ax = plt.subplots()

for backend, group in stats.groupby("backend"):
    group = group.sort_values("n_restarts")
    (line,) = ax.plot(group["n_restarts"], group["mean"], marker="o", label=backend)

ax.set_xlabel("# restarts")
ax.set_ylabel("mean run time (s)")
ax.legend()
plt.show()
