# %%
#
# This benchmark compares the speed of PCA solvers on datasets of different
# sizes in order to determine the best solver to select by default via the
# "auto" heuristic.
#
# Note: we do not control for the accuracy of the solvers: we assume that all
# solvers yield transformed data with similar explained variance. This
# assumption is generally true, except for the randomized solver that might
# require more power iterations.
#
# We generate synthetic data with dimensions that are useful to plot:
# - time vs n_samples for a fixed n_features and,
# - time vs n_features for a fixed n_samples for a fixed n_features.
import itertools
from math import log10
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn import config_context
from sklearn.decomposition import PCA

REF_DIMS = [100, 1000, 10_000]
data_shapes = []
for ref_dim in REF_DIMS:
    data_shapes.extend([(ref_dim, 10**i) for i in range(1, 8 - int(log10(ref_dim)))])
    data_shapes.extend(
        [(ref_dim, 3 * 10**i) for i in range(1, 8 - int(log10(ref_dim)))]
    )
    data_shapes.extend([(10**i, ref_dim) for i in range(1, 8 - int(log10(ref_dim)))])
    data_shapes.extend(
        [(3 * 10**i, ref_dim) for i in range(1, 8 - int(log10(ref_dim)))]
    )

# Remove duplicates:
data_shapes = sorted(set(data_shapes))

print("Generating test datasets...")
rng = np.random.default_rng(0)
datasets = [rng.normal(size=shape) for shape in data_shapes]


# %%
def measure_one(data, n_components, solver, method_name="fit"):
    print(
        f"Benchmarking {solver=!r}, {n_components=}, {method_name=!r} on data with"
        f" shape {data.shape}"
    )
    pca = PCA(n_components=n_components, svd_solver=solver, random_state=0)
    timings = []
    elapsed = 0
    method = getattr(pca, method_name)
    with config_context(assume_finite=True):
        while elapsed < 0.5:
            tic = perf_counter()
            method(data)
            duration = perf_counter() - tic
            timings.append(duration)
            elapsed += duration
    return np.median(timings)


SOLVERS = ["full", "covariance_eigh", "arpack", "randomized", "auto"]
measurements = []
for data, n_components, method_name in itertools.product(
    datasets, [2, 50], ["fit", "fit_transform"]
):
    if n_components >= min(data.shape):
        continue
    for solver in SOLVERS:
        if solver == "covariance_eigh" and data.shape[1] > 5000:
            # Too much memory and too slow.
            continue
        if solver in ["arpack", "full"] and log10(data.size) > 7:
            # Too slow, in particular for the full solver.
            continue
        time = measure_one(data, n_components, solver, method_name=method_name)
        measurements.append(
            {
                "n_components": n_components,
                "n_samples": data.shape[0],
                "n_features": data.shape[1],
                "time": time,
                "solver": solver,
                "method_name": method_name,
            }
        )
measurements = pd.DataFrame(measurements)
measurements.to_csv("bench_pca_solvers.csv", index=False)

# %%
all_method_names = measurements["method_name"].unique()
all_n_components = measurements["n_components"].unique()

for method_name in all_method_names:
    fig, axes = plt.subplots(
        figsize=(16, 16),
        nrows=len(REF_DIMS),
        ncols=len(all_n_components),
        sharey=True,
        constrained_layout=True,
    )
    fig.suptitle(f"Benchmarks for PCA.{method_name}, varying n_samples", fontsize=16)

    for row_idx, ref_dim in enumerate(REF_DIMS):
        for n_components, ax in zip(all_n_components, axes[row_idx]):
            for solver in SOLVERS:
                if solver == "auto":
                    style_kwargs = dict(linewidth=2, color="black", style="--")
                else:
                    style_kwargs = dict(style="o-")
                ax.set(
                    title=f"n_components={n_components}, n_features={ref_dim}",
                    ylabel="time (s)",
                )
                measurements.query(
                    "n_components == @n_components and n_features == @ref_dim"
                    " and solver == @solver and method_name == @method_name"
                ).plot.line(
                    x="n_samples",
                    y="time",
                    label=solver,
                    logx=True,
                    logy=True,
                    ax=ax,
                    **style_kwargs,
                )
# %%
for method_name in all_method_names:
    fig, axes = plt.subplots(
        figsize=(16, 16),
        nrows=len(REF_DIMS),
        ncols=len(all_n_components),
        sharey=True,
    )
    fig.suptitle(f"Benchmarks for PCA.{method_name}, varying n_features", fontsize=16)

    for row_idx, ref_dim in enumerate(REF_DIMS):
        for n_components, ax in zip(all_n_components, axes[row_idx]):
            for solver in SOLVERS:
                if solver == "auto":
                    style_kwargs = dict(linewidth=2, color="black", style="--")
                else:
                    style_kwargs = dict(style="o-")
                ax.set(
                    title=f"n_components={n_components}, n_samples={ref_dim}",
                    ylabel="time (s)",
                )
                measurements.query(
                    "n_components == @n_components and n_samples == @ref_dim "
                    " and solver == @solver and method_name == @method_name"
                ).plot.line(
                    x="n_features",
                    y="time",
                    label=solver,
                    logx=True,
                    logy=True,
                    ax=ax,
                    **style_kwargs,
                )

# %%
