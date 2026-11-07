"""
=================================================
GPU-Accelerated Anomaly Detection with GMM
=================================================

Demonstrates GPU-accelerated anomaly detection using Gaussian Mixture Models
on the credit card fraud dataset (284K transactions, 29 features). A pipeline
of StandardScaler and GaussianMixture is tuned via GridSearchCV over the
number of mixture components. Per-sample log-likelihoods from the best model
are used as anomaly scores to separate fraudulent from legitimate
transactions.
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

from joblib import parallel_backend

import sklearn
from sklearn.datasets import fetch_openml
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load the dataset
# ----------------
# The credit card fraud dataset contains 284K transactions with 29
# numerical features (PCA-transformed for privacy). Only 492 transactions
# (0.17%) are fraudulent.

cc = fetch_openml(data_id=1597, as_frame=False)
X_np = cc.data.astype(np.float32)
y_labels = cc.target

is_fraud = y_labels == "1"
print(f"Samples: {X_np.shape[0]}, Features: {X_np.shape[1]}")
print(f"Legitimate: {(~is_fraud).sum()}, Fraud: {is_fraud.sum()}")

# %%
# Define the pipeline and parameter grid
# ---------------------------------------
# GaussianMixture requires ``init_params="random_from_data"`` for Array
# API compatibility. We search over the number of mixture components with
# ``reg_covar=1e-4`` for numerical stability.

def make_gmm_pipeline():
    return make_pipeline(
        StandardScaler(),
        GaussianMixture(
            init_params="random_from_data",
            reg_covar=1e-4,
            random_state=42,
        ),
    )

param_grid = {"gaussianmixture__n_components": [2, 5, 10, 15, 20, 25, 30]}
search = GridSearchCV(make_gmm_pipeline(), param_grid, cv=3)

# %%
# CPU / NumPy baseline
# --------------------

t0 = time.perf_counter()
search.fit(X_np)
elapsed_cpu = time.perf_counter() - t0

best_n_cpu = search.best_params_["gaussianmixture__n_components"]
print(f"CPU time:                {elapsed_cpu:.1f}s")
print(f"Best n_components (CPU): {best_n_cpu}")

scores_cpu = search.best_estimator_.score_samples(X_np)

# %%
# GPU / PyTorch run
# -----------------

X_torch = torch.asarray(X_np, device=device, dtype=torch.float32)

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

search_gpu = GridSearchCV(make_gmm_pipeline(), param_grid, cv=3)

t0 = time.perf_counter()
search_gpu.fit(X_torch)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0

best_n_gpu = search_gpu.best_params_["gaussianmixture__n_components"]
print(f"GPU time:                {elapsed_gpu:.1f}s")
print(f"Best n_components (GPU): {best_n_gpu}")

scores_gpu = search_gpu.best_estimator_.score_samples(X_torch)
scores_gpu_np = (
    np.asarray(scores_gpu.cpu()) if hasattr(scores_gpu, "cpu") else scores_gpu
)

sklearn.set_config(array_api_dispatch=False)

# %%
# Compare timings
# ---------------

speedup = elapsed_cpu / elapsed_gpu if elapsed_gpu > 0 else float("inf")

fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(["NumPy / CPU", f"PyTorch / {device.upper()}"],
              [elapsed_cpu, elapsed_gpu])
ax.set_ylabel("Wall-clock time (seconds)")
ax.set_title(
    f"GMM GridSearchCV on Credit Card Fraud\nSpeedup: {speedup:.1f}x"
)
ax.bar_label(bars, fmt="%.1fs")
plt.tight_layout()
plt.show()

# %%
# Anomaly score distributions
# ---------------------------
# Low log-likelihood indicates an unusual transaction. The histogram shows
# that fraudulent transactions tend to concentrate in the low-likelihood
# tail.

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(scores_gpu_np[~is_fraud], bins=100, alpha=0.6, label="Legitimate", density=True)
ax.hist(scores_gpu_np[is_fraud], bins=100, alpha=0.6, label="Fraud", density=True)
ax.set_xlabel("Log-likelihood (anomaly score)")
ax.set_ylabel("Density")
ax.set_title("GMM Anomaly Score Distribution")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# ROC-AUC of the anomaly detector
# --------------------------------
# Using negative log-likelihood as the anomaly score (higher = more
# anomalous), we compute the ROC-AUC to quantify separability.

y_binary = is_fraud.astype(int)
auc = roc_auc_score(y_binary, -scores_gpu_np)
print(f"ROC-AUC (anomaly detection): {auc:.4f}")
