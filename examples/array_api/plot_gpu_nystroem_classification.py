"""
=====================================================
Non-Linear Classification at Scale with Nystroem
=====================================================

.. note::

    This example does not currently demonstrate a meaningful GPU speedup.
    The pipeline's compute is not large enough to overcome GPU dispatch
    overhead. This example needs to be reworked or replaced.

Demonstrates kernel approximation on the Covertype dataset (581K samples).
A Nystroem RBF feature map followed by RidgeClassifierCV inside
CalibratedClassifierCV delivers non-linear classification with calibrated
probabilities.

This contrasts with the linear PCA + LogisticRegression approach shown in
``plot_gpu_ridge_tuning.py`` on the same dataset: kernel features capture
non-linear decision boundaries that a linear model cannot.
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
from sklearn.calibration import CalibratedClassifierCV
from sklearn.datasets import fetch_covtype
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import RidgeClassifierCV
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load the dataset
# ----------------
# Covertype has 581K samples and 54 features -- large enough for the
# Nystroem kernel matrix computation to be a genuine bottleneck on CPU.

X_np, y_np = fetch_covtype(return_X_y=True)

X_train_np, X_test_np, y_train_np, y_test_np = train_test_split(
    X_np, y_np, test_size=0.2, random_state=42, stratify=y_np
)
print(f"Train: {X_train_np.shape[0]}, Test: {X_test_np.shape[0]}")

# %%
# Define the pipeline
# -------------------
# StandardScaler normalises features. Nystroem constructs an RBF kernel
# feature map with 1000 components, projecting each sample into a
# 1000-dimensional space that approximates the full RBF kernel matrix.
# CalibratedClassifierCV wraps RidgeClassifierCV to produce calibrated
# probabilities via temperature scaling (the only calibration method
# compatible with Array API).

pipe = make_pipeline(
    StandardScaler(),
    Nystroem(kernel="rbf", n_components=1000, random_state=42),
    CalibratedClassifierCV(
        RidgeClassifierCV(alphas=np.logspace(-3, 3, 7)),
        method="temperature",
    ),
)

# %%
# CPU / NumPy baseline
# --------------------

t0 = time.perf_counter()
pipe.fit(X_train_np, y_train_np)
y_pred_cpu = pipe.predict(X_test_np)
elapsed_cpu = time.perf_counter() - t0

acc_cpu = accuracy_score(y_test_np, y_pred_cpu)
print(f"CPU time:     {elapsed_cpu:.1f}s")
print(f"CPU accuracy: {acc_cpu:.4f}")

# %%
# GPU / PyTorch run
# -----------------

X_train_torch = torch.asarray(X_train_np, device=device, dtype=torch.float32)
y_train_torch = torch.asarray(y_train_np, device=device, dtype=torch.float32)
X_test_torch = torch.asarray(X_test_np, device=device, dtype=torch.float32)
y_test_torch = torch.asarray(y_test_np, device=device, dtype=torch.float32)

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

pipe_gpu = make_pipeline(
    StandardScaler(),
    Nystroem(kernel="rbf", n_components=1000, random_state=42),
    CalibratedClassifierCV(
        RidgeClassifierCV(alphas=np.logspace(-3, 3, 7)),
        method="temperature",
    ),
)

t0 = time.perf_counter()
pipe_gpu.fit(X_train_torch, y_train_torch)
y_pred_gpu = pipe_gpu.predict(X_test_torch)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0

acc_gpu = float(accuracy_score(y_test_torch, y_pred_gpu))
print(f"GPU time:     {elapsed_gpu:.1f}s")
print(f"GPU accuracy: {acc_gpu:.4f}")

sklearn.set_config(array_api_dispatch=False)

# %%
# Results
# -------
# Both backends should produce comparable accuracy (small numerical
# differences are expected from float32 arithmetic).

speedup = elapsed_cpu / elapsed_gpu if elapsed_gpu > 0 else float("inf")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

bars = axes[0].bar(["NumPy / CPU", f"PyTorch / {device.upper()}"],
                   [elapsed_cpu, elapsed_gpu])
axes[0].set_ylabel("Wall-clock time (seconds)")
axes[0].set_title(
    f"Nystroem + RidgeClassifierCV on Covertype\nSpeedup: {speedup:.1f}x"
)
axes[0].bar_label(bars, fmt="%.1fs")

axes[1].bar(["CPU", "GPU"], [acc_cpu, acc_gpu])
axes[1].set_ylabel("Accuracy")
axes[1].set_title("Test Accuracy")
axes[1].set_ylim(0.8, 1.0)

plt.tight_layout()
plt.show()

# %%
# Classification report (GPU model)
# ----------------------------------

y_pred_np = np.asarray(y_pred_gpu.cpu()) if hasattr(y_pred_gpu, "cpu") else y_pred_gpu
y_test_report = np.asarray(
    y_test_torch.cpu()) if hasattr(y_test_torch, "cpu") else y_test_torch
print(classification_report(y_test_report.astype(int), y_pred_np.astype(int)))
