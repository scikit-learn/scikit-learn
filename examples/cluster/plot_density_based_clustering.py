"""
==============================================
Density-based clustering: DBSCAN, OPTICS, HDBSCAN
==============================================

This example compares three density-based clustering algorithms on the same
synthetic dataset: DBSCAN, OPTICS, and HDBSCAN. All three algorithms discover
clusters based on local point density and can label sparse points as noise.

DBSCAN uses a single neighborhood radius (``eps``), so it can struggle when
cluster densities vary. OPTICS orders points by reachability distance so that
one fit can be used to extract clusters at multiple ``eps`` values. HDBSCAN
builds a hierarchy across density levels and extracts stable clusters without
requiring a global ``eps``.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate data
# -------------
import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import DBSCAN, HDBSCAN, OPTICS, cluster_optics_dbscan
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler

centers = [[-0.85, -0.85], [-0.85, 0.85], [3, 3], [3, -3]]
X, labels_true = make_blobs(
    n_samples=750,
    centers=centers,
    cluster_std=[0.2, 0.35, 1.35, 1.35],
    random_state=0,
)
X = StandardScaler(copy=True).fit_transform(X)


# %%
# Helper to plot cluster labels
# -----------------------------


def plot_labels(ax, X, labels, title, probabilities=None):
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    n_noise = int(np.sum(labels == -1))

    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]
        mask = labels == k
        marker = "o" if k != -1 else "x"
        if probabilities is None:
            sizes = 14
        else:
            sizes = 5 + 40 * probabilities[mask]
            if k == -1:
                sizes = 14

        ax.scatter(
            X[mask, 0],
            X[mask, 1],
            c=[col],
            s=sizes,
            marker=marker,
            alpha=0.8,
        )

    ax.set_title(f"{title}\n{n_clusters} clusters, {n_noise} noise points")


# %%
# Ground truth
# ------------
fig, ax = plt.subplots(figsize=(8, 5))
unique_labels = set(labels_true)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    mask = labels_true == k
    ax.scatter(X[mask, 0], X[mask, 1], c=[col], s=14, label=f"Cluster {k}")
ax.set_title("Ground truth clusters")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# DBSCAN: sensitivity to a single eps
# -----------------------------------
# Two different eps values show how a single global radius can under-cluster
# sparse regions or merge tight clusters.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# eps tuned for tight clusters - sparse clusters fragment into noise
small_eps = 0.15
db_small = DBSCAN(eps=small_eps, min_samples=5).fit(X)
plot_labels(ax1, X, db_small.labels_, f"DBSCAN eps={small_eps:.2f}")

# eps tuned for sparse clusters - tight clusters merge incorrectly
large_eps = 0.60
db_large = DBSCAN(eps=large_eps, min_samples=5).fit(X)
plot_labels(ax2, X, db_large.labels_, f"DBSCAN eps={large_eps:.2f}")

fig.suptitle("DBSCAN: no single eps fits all densities", fontsize=13)
plt.tight_layout()
plt.show()

# %%
# OPTICS: reachability plot
# -------------------------
optics = OPTICS(min_samples=10, xi=0.05, min_cluster_size=0.05)
optics.fit(X)

reachability = optics.reachability_[optics.ordering_]
labels_optics = optics.labels_[optics.ordering_]
space = np.arange(len(X))

fig, ax = plt.subplots(figsize=(10, 4))
unique_labels = set(labels_optics)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    mask = labels_optics == k
    if k == -1:
        col = [0, 0, 0, 0.4]
    ax.plot(space[mask], reachability[mask], ".", color=col, markersize=2)

ax.axhline(0.20, linestyle="--", color="k", linewidth=1, label="eps=0.20")
ax.axhline(0.75, linestyle="--", color="r", linewidth=1, label="eps=0.75")
ax.set_ylabel("Reachability distance")
ax.set_xlabel("Point ordering")
ax.set_title("OPTICS reachability plot")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# OPTICS: cluster extraction at two thresholds
# --------------------------------------------
labels_020 = cluster_optics_dbscan(
    reachability=optics.reachability_,
    core_distances=optics.core_distances_,
    ordering=optics.ordering_,
    eps=0.20,
)
labels_075 = cluster_optics_dbscan(
    reachability=optics.reachability_,
    core_distances=optics.core_distances_,
    ordering=optics.ordering_,
    eps=0.75,
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plot_labels(ax1, X, labels_020, "OPTICS eps=0.20")
plot_labels(ax2, X, labels_075, "OPTICS eps=0.75")
fig.suptitle("OPTICS: same fit, different eps thresholds", fontsize=13)
plt.tight_layout()
plt.show()

# %%
# HDBSCAN: automatic multi-scale clustering
# -----------------------------------------
# HDBSCAN does not require a global eps. Marker size reflects membership
# probability, so low-confidence points appear smaller.

hdb = HDBSCAN(min_cluster_size=15, min_samples=5, copy=True)
hdb.fit(X)

fig, ax = plt.subplots(figsize=(8, 5))
plot_labels(
    ax,
    X,
    hdb.labels_,
    "HDBSCAN min_cluster_size=15\nMarker size = membership probability",
    probabilities=hdb.probabilities_,
)
plt.tight_layout()
plt.show()

# %%
# Side-by-side comparison
# -----------------------
# DBSCAN requires manual tuning, OPTICS helps choose a threshold, and HDBSCAN
# extracts clusters automatically from the density hierarchy.

best_eps = 0.30
db_best = DBSCAN(eps=best_eps, min_samples=5).fit(X)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
plot_labels(axes[0], X, db_best.labels_, f"DBSCAN eps={best_eps:.2f}")
plot_labels(axes[1], X, labels_075, "OPTICS eps=0.75")
plot_labels(axes[2], X, hdb.labels_, "HDBSCAN min_cluster_size=15")
fig.suptitle("Comparing DBSCAN, OPTICS, and HDBSCAN on varying densities", fontsize=12)
plt.tight_layout()
plt.show()
