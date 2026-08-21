"""
===============================================
Semantic Topic Discovery via Embedding Clustering
===============================================

Demonstrates GPU-accelerated unsupervised topic discovery on a text corpus.
Sentence embeddings from ``all-MiniLM-L6-v2`` are reduced with PCA and
clustered with GaussianMixture, with the number of topics selected via
GridSearchCV. The entire pipeline -- from dimensionality reduction through
EM iterations to cluster evaluation -- runs on GPU when Array API dispatch
is enabled.

This example requires ``sentence-transformers`` (``pip install
sentence-transformers``).
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

try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    print(
        "This example requires sentence-transformers."
        " Install it with: pip install sentence-transformers"
    )
    sys.exit(0)

import sklearn
from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.cluster import calinski_harabasz_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

# %%
# Load dataset and compute embeddings
# ------------------------------------

newsgroups = fetch_20newsgroups(subset="all")
all_texts = newsgroups.data
all_targets = newsgroups.target
target_names = newsgroups.target_names

print(f"Documents: {len(all_texts)}, True categories: {len(target_names)}")

encoder = SentenceTransformer("all-MiniLM-L6-v2", device=device)

embeddings_gpu = encoder.encode(
    all_texts, convert_to_tensor=True, batch_size=128, show_progress_bar=True
)
embeddings_np = np.asarray(
    embeddings_gpu.cpu() if hasattr(embeddings_gpu, "cpu") else embeddings_gpu
)
print(f"Embedding shape: {tuple(embeddings_gpu.shape)}")

# %%
# Define the clustering pipeline
# -------------------------------
# StandardScaler + PCA(50) + GaussianMixture with ``init_params="random_from_data"``
# (required for Array API). GridSearchCV selects the number of mixture
# components.

param_grid = {"gaussianmixture__n_components": [10, 15, 20, 25, 30]}

def make_topic_pipeline():
    return make_pipeline(
        StandardScaler(),
        PCA(svd_solver="full", n_components=50),
        GaussianMixture(
            init_params="random_from_data",
            random_state=42,
        ),
    )

# %%
# CPU / NumPy baseline
# --------------------

search_cpu = GridSearchCV(make_topic_pipeline(), param_grid, cv=3)

t0 = time.perf_counter()
search_cpu.fit(embeddings_np)
elapsed_cpu = time.perf_counter() - t0

labels_cpu = search_cpu.best_estimator_.predict(embeddings_np)
best_n_cpu = search_cpu.best_params_["gaussianmixture__n_components"]
print(f"CPU time:                {elapsed_cpu:.1f}s")
print(f"Best n_components (CPU): {best_n_cpu}")

# %%
# GPU / PyTorch run
# -----------------

if device == "cuda":
    torch.zeros(1, device="cuda")

sklearn.set_config(array_api_dispatch=True)

search_gpu = GridSearchCV(make_topic_pipeline(), param_grid, cv=3)

t0 = time.perf_counter()
search_gpu.fit(embeddings_gpu)
if device == "cuda":
    torch.cuda.synchronize()
elapsed_gpu = time.perf_counter() - t0

labels_gpu = search_gpu.best_estimator_.predict(embeddings_gpu)
labels_gpu_np = np.asarray(
    labels_gpu.cpu() if hasattr(labels_gpu, "cpu") else labels_gpu
)
best_n_gpu = search_gpu.best_params_["gaussianmixture__n_components"]
print(f"GPU time:                {elapsed_gpu:.1f}s")
print(f"Best n_components (GPU): {best_n_gpu}")

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
    f"GMM Topic Discovery on 20 Newsgroups\nSpeedup: {speedup:.1f}x"
)
ax.bar_label(bars, fmt="%.1fs")
plt.tight_layout()
plt.show()

# %%
# Evaluate discovered topics against true categories
# ---------------------------------------------------
# Since the 20 Newsgroups dataset has ground-truth labels, we can measure
# how well the discovered clusters align with them.

ari = adjusted_rand_score(all_targets, labels_gpu_np)
nmi = normalized_mutual_info_score(all_targets, labels_gpu_np)
ch = calinski_harabasz_score(embeddings_np, labels_gpu_np)

print(f"Adjusted Rand Index:           {ari:.4f}")
print(f"Normalized Mutual Information: {nmi:.4f}")
print(f"Calinski-Harabasz Score:       {ch:.1f}")

# %%
# 2D PCA visualisation of discovered topics
# ------------------------------------------
# Project embeddings to 2D with PCA for a scatter plot coloured by the
# discovered cluster assignments.

pca_2d = PCA(n_components=2, svd_solver="full")
coords = pca_2d.fit_transform(embeddings_np)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

scatter0 = axes[0].scatter(
    coords[:, 0], coords[:, 1],
    c=labels_gpu_np, cmap="tab20", s=1, alpha=0.5,
)
axes[0].set_title(f"Discovered Topics ({best_n_gpu} clusters)")
axes[0].set_xlabel("PC 1")
axes[0].set_ylabel("PC 2")
plt.colorbar(scatter0, ax=axes[0], label="Cluster")

scatter1 = axes[1].scatter(
    coords[:, 0], coords[:, 1],
    c=all_targets, cmap="tab20", s=1, alpha=0.5,
)
axes[1].set_title("True Newsgroup Categories")
axes[1].set_xlabel("PC 1")
axes[1].set_ylabel("PC 2")
plt.colorbar(scatter1, ax=axes[1], label="Category")

plt.tight_layout()
plt.show()

# %%
# Topic composition
# -----------------
# For each discovered cluster, show the dominant true newsgroup category
# to give a sense of what "topic" each cluster captures.

print(f"\n{'Cluster':>8}  {'Size':>6}  Dominant Category")
print("-" * 50)
for cluster_id in range(best_n_gpu):
    mask = labels_gpu_np == cluster_id
    if mask.sum() == 0:
        continue
    category_counts = np.bincount(all_targets[mask], minlength=len(target_names))
    dominant = target_names[np.argmax(category_counts)]
    purity = category_counts.max() / mask.sum()
    print(f"{cluster_id:>8}  {mask.sum():>6}  {dominant} ({purity:.0%})")
