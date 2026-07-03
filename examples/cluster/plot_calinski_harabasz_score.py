"""
========================================================================
Comparing Calinski-Harabasz and Silhouette scores for cluster evaluation
========================================================================

This example compares two unsupervised cluster evaluation metrics —
:func:`~sklearn.metrics.calinski_harabasz_score` and
:func:`~sklearn.metrics.silhouette_score` — to illustrate when each is
appropriate for selecting the optimal number of clusters.

**Calinski-Harabasz score** (also known as the Variance Ratio Criterion)
is the ratio of the between-cluster dispersion to the within-cluster
dispersion. A higher score indicates better-defined clusters. The score is
unbounded above (no maximum value) and scales with the number of samples,
making direct cross-dataset comparisons unreliable — but it is
substantially faster to compute than the silhouette score, which matters
on large datasets.

**Silhouette score** measures how similar a sample is to its own cluster
compared to other clusters. It ranges from -1 (incorrect clustering) to
+1 (highly dense clustering), with values near 0 indicating overlapping
clusters.

The example is divided into three sections:

- **Well-separated blobs**: both metrics agree on the optimal
  ``n_clusters``.
- **Noisy / overlapping data**: the metrics can diverge; understanding
  why helps you choose the right one.
- **Runtime comparison**: Calinski-Harabasz is much faster on large
  datasets.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Well-separated clusters: both metrics agree
# -------------------------------------------
#
# We generate isotropic Gaussian blobs with a clear ground truth of 4
# clusters and sweep ``n_clusters`` from 2 to 8. Both metrics should peak
# at ``n_clusters=4``.

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.metrics import calinski_harabasz_score, silhouette_score

X_blobs, _ = make_blobs(n_samples=500, centers=4, cluster_std=0.8, random_state=42)

n_clusters_range = range(2, 9)
ch_scores = []
sil_scores = []

for k in n_clusters_range:
    labels = KMeans(n_clusters=k, random_state=42, n_init="auto").fit_predict(X_blobs)
    ch_scores.append(calinski_harabasz_score(X_blobs, labels))
    sil_scores.append(silhouette_score(X_blobs, labels))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(list(n_clusters_range), ch_scores, marker="o", color="steelblue")
axes[0].set_xlabel("Number of clusters")
axes[0].set_ylabel("Calinski-Harabasz score")
axes[0].set_title("Calinski-Harabasz (higher is better, unbounded)")
axes[0].axvline(x=4, color="red", linestyle="--", label="True k=4")
axes[0].legend()

axes[1].plot(list(n_clusters_range), sil_scores, marker="o", color="darkorange")
axes[1].set_xlabel("Number of clusters")
axes[1].set_ylabel("Silhouette score")
axes[1].set_title("Silhouette (higher is better, range [-1, 1])")
axes[1].axvline(x=4, color="red", linestyle="--", label="True k=4")
axes[1].legend()

fig.suptitle("Well-separated blobs: both metrics peak at the true k=4", fontsize=13)
plt.tight_layout()
plt.show()

# %%
# The Calinski-Harabasz score peaks sharply at ``k=4`` and then drops.
# The silhouette score behaves similarly. When clusters are well-separated
# and roughly convex, both metrics reliably identify the correct number of
# clusters.

# %%
# Noisy / overlapping data: metrics can disagree
# -----------------------------------------------
#
# When clusters overlap or data contains noise, the two metrics can give
# different signals. Here we add noise points and use tighter, overlapping
# blobs. Calinski-Harabasz is sensitive to cluster size imbalance and
# may prefer solutions with more, smaller groups because that reduces
# within-cluster variance. Silhouette penalises points that are close to
# the boundary of a neighbouring cluster, which can make it more
# conservative.

from sklearn.datasets import make_moons

X_moons, _ = make_moons(n_samples=300, noise=0.1, random_state=42)

# Add a small amount of uniform noise
rng = np.random.default_rng(0)
X_noisy = np.vstack([X_moons, rng.uniform(-1.5, 2.5, size=(30, 2))])

ch_noisy = []
sil_noisy = []

for k in n_clusters_range:
    labels = KMeans(n_clusters=k, random_state=42, n_init="auto").fit_predict(X_noisy)
    ch_noisy.append(calinski_harabasz_score(X_noisy, labels))
    sil_noisy.append(silhouette_score(X_noisy, labels))

fig, axes = plt.subplots(1, 3, figsize=(16, 4))

axes[0].scatter(X_noisy[:, 0], X_noisy[:, 1], s=10, alpha=0.6)
axes[0].set_title("Noisy two-moons dataset")

axes[1].plot(list(n_clusters_range), ch_noisy, marker="o", color="steelblue")
axes[1].set_xlabel("Number of clusters")
axes[1].set_ylabel("Calinski-Harabasz score")
axes[1].set_title("Calinski-Harabasz")

axes[2].plot(list(n_clusters_range), sil_noisy, marker="o", color="darkorange")
axes[2].set_xlabel("Number of clusters")
axes[2].set_ylabel("Silhouette score")
axes[2].set_title("Silhouette")

fig.suptitle(
    "Overlapping / non-convex data: metrics may suggest different k", fontsize=13
)
plt.tight_layout()
plt.show()

# %%
# With non-convex geometry (moons) and added noise, the two metrics may
# peak at different values of ``k``. Neither is "wrong" — they capture
# different aspects of cluster structure. In such cases, inspecting the
# actual cluster assignments alongside both scores is recommended.

# %%
# Runtime comparison on large data
# ---------------------------------
#
# Calinski-Harabasz is :math:`O(n \cdot d)` (a single pass over the
# data), while Silhouette is :math:`O(n^2)` (pairwise distances). On
# large datasets the difference is significant.

import time

from sklearn.datasets import make_blobs

X_large, _ = make_blobs(n_samples=5_000, centers=5, random_state=42)
labels_large = KMeans(n_clusters=5, random_state=42, n_init="auto").fit_predict(
    X_large
)

t0 = time.perf_counter()
for _ in range(10):
    calinski_harabasz_score(X_large, labels_large)
ch_time = (time.perf_counter() - t0) / 10

t0 = time.perf_counter()
for _ in range(10):
    silhouette_score(X_large, labels_large, sample_size=1_000, random_state=42)
sil_time = (time.perf_counter() - t0) / 10

fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(
    ["Calinski-Harabasz", "Silhouette\n(sample_size=1000)"],
    [ch_time * 1000, sil_time * 1000],
    color=["steelblue", "darkorange"],
)
ax.bar_label(bars, fmt="{:.1f} ms")
ax.set_ylabel("Time per call (ms)")
ax.set_title(f"Runtime on {len(X_large):,} samples (averaged over 10 runs)")
plt.tight_layout()
plt.show()

# %%
# Even when silhouette is approximated with ``sample_size=1000``, it is
# noticeably slower than the exact Calinski-Harabasz computation. For
# datasets with tens of thousands of samples, prefer Calinski-Harabasz
# for fast hyper-parameter sweeps and use silhouette as a secondary check
# on candidate solutions.
#
# .. rubric:: Summary
#
# +---------------------------+-------------------------------+---------------------------+
# | Property                  | Calinski-Harabasz             | Silhouette                |
# +===========================+===============================+===========================+
# | Range                     | [0, ∞) — higher is better    | [-1, 1] — higher is better|
# +---------------------------+-------------------------------+---------------------------+
# | Complexity                | O(n · d) — fast               | O(n²) — slow              |
# +---------------------------+-------------------------------+---------------------------+
# | Sensitive to cluster size | Yes (prefers balanced sizes)  | Less so                   |
# +---------------------------+-------------------------------+---------------------------+
# | Works with non-convex     | Poor                          | Poor (density-based       |
# | clusters                  |                               | metrics preferred)        |
# +---------------------------+-------------------------------+---------------------------+
# | Good for                  | Fast sweeps, large data       | Interpretable [-1,1] score|
# +---------------------------+-------------------------------+---------------------------+
