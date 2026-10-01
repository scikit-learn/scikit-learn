"""
=============================================
GPU-Accelerated Hyperparameter Search
=============================================

Demonstrates how scikit-learn's Array API support accelerates a full
hyperparameter search on a large tabular dataset. A pipeline of
StandardScaler, PCA, and LogisticRegression is tuned with GridSearchCV over
25 parameter combinations on the Covertype dataset (581K samples). When run
on a CUDA GPU, all linear algebra operations (SVD for PCA, L-BFGS iterations
for logistic regression) are offloaded to the GPU, yielding significant
wall-clock speedups with zero code changes to the pipeline definition.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
import time

os.environ["SCIPY_ARRAY_API"] = "1"

import matplotlib.pyplot as plt
import numpy as np

try:
    import torch
except ImportError:
    print("This example requires PyTorch. Install it with: pip install torch")
    sys.exit(0)

import sklearn
from sklearn.datasets import fetch_covtype
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load the dataset
# ----------------
# The Covertype dataset has 581K samples and 54 features, making it large
# enough for GPU acceleration to matter.

X_np, y_np = fetch_covtype(return_X_y=True)
print(f"Dataset: {X_np.shape[0]} samples, {X_np.shape[1]} features")

# %%
# Define the pipeline and parameter grid
# ---------------------------------------
# PCA requires ``svd_solver="full"`` and LogisticRegression requires
# ``solver="lbfgs"`` for Array API compatibility.

pipe = make_pipeline(
    StandardScaler(),
    PCA(svd_solver="full"),
    LogisticRegression(solver="lbfgs", max_iter=1000),
)

param_grid = {
    "pca__n_components": [10, 20, 30, 40, 50],
    "logisticregression__C": np.logspace(-2, 2, 5).tolist(),
}

search = GridSearchCV(pipe, param_grid, cv=3, n_jobs=32)

# %%
# CPU / NumPy baseline
# --------------------
# Fit the grid search using standard NumPy arrays on the CPU.

t0 = time.perf_counter()
search.fit(X_np, y_np)
elapsed_cpu = time.perf_counter() - t0

best_params_cpu = search.best_params_
best_score_cpu = search.best_score_
print(f"CPU time:       {elapsed_cpu:.1f}s")
print(f"Best CV score:  {best_score_cpu:.4f}")
print(f"Best params:    {best_params_cpu}")

# %%
# GPU / PyTorch run
# -----------------
# Convert data to PyTorch tensors on the target device, enable Array API
# dispatch, and re-run the same grid search. The pipeline definition is
# unchanged.

X_torch = torch.asarray(X_np, device=device, dtype=torch.float32)
y_torch = torch.asarray(y_np, device=device, dtype=torch.float32)

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

search_gpu = GridSearchCV(pipe, param_grid, cv=3)

t0 = time.perf_counter()
search_gpu.fit(X_torch, y_torch)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0

best_params_gpu = search_gpu.best_params_
best_score_gpu = search_gpu.best_score_
print(f"GPU time:       {elapsed_gpu:.1f}s")
print(f"Best CV score:  {best_score_gpu:.4f}")
print(f"Best params:    {best_params_gpu}")

sklearn.set_config(array_api_dispatch=False)

# %%
# Compare timings
# ---------------
# The speedup comes from GPU-accelerated BLAS operations inside PCA (SVD)
# and LogisticRegression (L-BFGS), multiplied across 75 pipeline fits
# (25 param combos x 3 CV folds).

speedup = elapsed_cpu / elapsed_gpu if elapsed_gpu > 0 else float("inf")

fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(["NumPy / CPU", "PyTorch / " + device.upper()],
              [elapsed_cpu, elapsed_gpu])
ax.set_ylabel("Wall-clock time (seconds)")
ax.set_title(f"GridSearchCV: 25 configs x 3 folds on Covertype\nSpeedup: {speedup:.1f}x")
ax.bar_label(bars, fmt="%.1fs")
plt.tight_layout()
plt.show()

