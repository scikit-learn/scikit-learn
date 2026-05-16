"""
=================================================================
Comparing DBSCAN, HDBSCAN, and OPTICS for density-based clustering
=================================================================

DBSCAN, HDBSCAN, and OPTICS are three density-based clustering algorithms
available in scikit-learn. They share a common idea -- finding regions of
high sample density and grouping them into clusters while labeling sparse
regions as noise -- but they differ in how they discover clusters:

- :class:`~sklearn.cluster.DBSCAN` uses a single density threshold (``eps``)
  to define neighborhoods. It is fast and simple, but it requires choosing
  a value of ``eps`` that fits the data.
- :class:`~sklearn.cluster.OPTICS` builds a *reachability plot* that encodes
  clusterings at every density level. Clusters can be extracted automatically
  via the Xi method or via a DBSCAN-like horizontal cut on the reachability.
- :class:`~sklearn.cluster.HDBSCAN` extends DBSCAN by considering all values
  of ``eps`` simultaneously and condensing the result into a hierarchy of
  clusters. It is scale-invariant and adapts to clusters of varying density.

This example compares them on two datasets to highlight where they agree
and where they diverge.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Uniform-density blobs: all three algorithms agree
# -------------------------------------------------
#
# We first generate three isotropic Gaussian blobs of roughly equal density
# using :func:`~sklearn.datasets.make_blobs`. When density is uniform and
# the hyperparameters are well tuned, the three algorithms produce visually
# equivalent partitions. The Adjusted Rand Index (ARI) confirms that each
# clustering closely matches the ground-truth labels.

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import DBSCAN, HDBSCAN, OPTICS, cluster_optics_dbscan
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(
    n_samples=750, centers=centers, cluster_std=0.4, random_state=0
)
X = StandardScaler().fit_transform(X)

estimators = {
    "DBSCAN": DBSCAN(eps=0.3, min_samples=10),
    "HDBSCAN": HDBSCAN(min_cluster_size=20),
    "OPTICS": OPTICS(min_samples=10, xi=0.05, min_cluster_size=0.05),
}

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, (name, est) in zip(axes, estimators.items()):
    labels = est.fit_predict(X)
    ari = adjusted_rand_score(labels_true, labels)
    n_noise = int(np.sum(labels == -1))
    ax.scatter(X[:, 0], X[:, 1], c=labels, cmap="Spectral", s=10)
    ax.set_title(f"{name}\nARI={ari:.2f}, noise={n_noise}")
    ax.set_xticks([])
    ax.set_yticks([])
plt.tight_layout()
plt.show()

# %%
# All three algorithms recover the three clusters with similar quality.
# Noisy samples are assigned the label ``-1``. When the data looks like
# this, DBSCAN is usually the right tool: it is the simplest and fastest
# of the three.

# %%
# Varying-density data: where DBSCAN struggles
# --------------------------------------------
#
# Density-based clustering becomes interesting when clusters have
# different intrinsic densities. We build a 6-cluster dataset where the
# Gaussian standard deviations vary by an order of magnitude, so a single
# value of ``eps`` cannot fit every cluster at once.

np.random.seed(0)
n_points = 250
C1 = [-5, -2] + 0.8 * np.random.randn(n_points, 2)
C2 = [4, -1] + 0.1 * np.random.randn(n_points, 2)
C3 = [1, -2] + 0.2 * np.random.randn(n_points, 2)
C4 = [-2, 3] + 0.3 * np.random.randn(n_points, 2)
C5 = [3, -2] + 1.6 * np.random.randn(n_points, 2)
C6 = [5, 6] + 2.0 * np.random.randn(n_points, 2)
X_var = np.vstack((C1, C2, C3, C4, C5, C6))

estimators_var = {
    "DBSCAN (eps=0.5)": DBSCAN(eps=0.5, min_samples=50),
    "OPTICS (Xi)": OPTICS(min_samples=50, xi=0.05, min_cluster_size=0.05),
    "HDBSCAN": HDBSCAN(min_cluster_size=50),
}

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, (name, est) in zip(axes, estimators_var.items()):
    labels = est.fit_predict(X_var)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int(np.sum(labels == -1))
    ax.scatter(X_var[:, 0], X_var[:, 1], c=labels, cmap="Spectral", s=8)
    ax.set_title(f"{name}\nclusters={n_clusters}, noise={n_noise}")
    ax.set_xticks([])
    ax.set_yticks([])
plt.tight_layout()
plt.show()

# %%
# With a single ``eps``, DBSCAN can either merge nearby dense clusters or
# fragment the sparse ones. OPTICS (with the Xi extraction method) and
# HDBSCAN both discover the multi-density structure without manual
# ``eps`` tuning.

# %%
# OPTICS reachability plot
# ------------------------
#
# A defining feature of OPTICS is the *reachability plot*: for each point,
# its reachability distance is plotted in the order OPTICS visits points.
# Valleys in this plot correspond to clusters. Cutting the plot
# horizontally at different heights recovers DBSCAN-like clusterings at
# different values of ``eps``, all from the same OPTICS fit.

optics = OPTICS(min_samples=50, xi=0.05, min_cluster_size=0.05).fit(X_var)
space = np.arange(len(X_var))
reachability = optics.reachability_[optics.ordering_]
labels_xi = optics.labels_[optics.ordering_]

labels_05 = cluster_optics_dbscan(
    reachability=optics.reachability_,
    core_distances=optics.core_distances_,
    ordering=optics.ordering_,
    eps=0.5,
)
labels_2 = cluster_optics_dbscan(
    reachability=optics.reachability_,
    core_distances=optics.core_distances_,
    ordering=optics.ordering_,
    eps=2.0,
)

plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(2, 2)
ax_reach = plt.subplot(gs[0, :])
ax_05 = plt.subplot(gs[1, 0])
ax_2 = plt.subplot(gs[1, 1])

n_xi_clusters = labels_xi.max() + 1
colors = plt.cm.Spectral(np.linspace(0, 1, max(n_xi_clusters, 1)))
for k in range(n_xi_clusters):
    mask = labels_xi == k
    ax_reach.plot(space[mask], reachability[mask], ".", color=colors[k], alpha=0.5)
ax_reach.plot(space[labels_xi == -1], reachability[labels_xi == -1], "k.", alpha=0.3)
ax_reach.axhline(0.5, color="k", ls="-", alpha=0.5, label="eps=0.5")
ax_reach.axhline(2.0, color="k", ls="--", alpha=0.5, label="eps=2.0")
ax_reach.set_ylabel("Reachability distance")
ax_reach.set_title("OPTICS reachability plot (sorted by traversal order)")
ax_reach.legend(loc="upper right")

ax_05.scatter(X_var[:, 0], X_var[:, 1], c=labels_05, cmap="Spectral", s=8)
ax_05.set_title("Horizontal cut at eps=0.5")
ax_05.set_xticks([])
ax_05.set_yticks([])

ax_2.scatter(X_var[:, 0], X_var[:, 1], c=labels_2, cmap="Spectral", s=8)
ax_2.set_title("Horizontal cut at eps=2.0")
ax_2.set_xticks([])
ax_2.set_yticks([])

plt.tight_layout()
plt.show()

# %%
# The reachability plot turns clustering into a one-dimensional problem:
# any horizontal cut yields a valid clustering. The Xi method extracts a
# canonical labelling automatically from the shape of the valleys, while
# the two cuts shown above recover the families of clusterings that
# DBSCAN would have produced with ``eps=0.5`` and ``eps=2.0``.

# %%
# HDBSCAN scale invariance
# ------------------------
#
# DBSCAN's ``eps`` is tied to the absolute scale of the data: rescaling
# the inputs by a constant factor changes which clusters appear. HDBSCAN
# instead extracts clusters from a hierarchy and is invariant to such
# rescaling.

centers = [[1, 1], [-1, -1], [1.5, -1.5]]
X_scale, _ = make_blobs(
    n_samples=750, centers=centers, cluster_std=[0.4, 0.1, 0.75], random_state=0
)

fig, axes = plt.subplots(2, 3, figsize=(12, 7))
scales = [0.5, 1.0, 3.0]
for col, scale in enumerate(scales):
    X_s = X_scale * scale

    dbs = DBSCAN(eps=0.3).fit(X_s)
    axes[0, col].scatter(X_s[:, 0], X_s[:, 1], c=dbs.labels_, cmap="Spectral", s=8)
    axes[0, col].set_title(f"DBSCAN (eps=0.3), scale={scale}")
    axes[0, col].set_xticks([])
    axes[0, col].set_yticks([])

    hdb = HDBSCAN().fit(X_s)
    axes[1, col].scatter(X_s[:, 0], X_s[:, 1], c=hdb.labels_, cmap="Spectral", s=8)
    axes[1, col].set_title(f"HDBSCAN, scale={scale}")
    axes[1, col].set_xticks([])
    axes[1, col].set_yticks([])

plt.tight_layout()
plt.show()

# %%
# A fixed ``eps`` works at one scale and breaks at the others, while
# HDBSCAN produces an essentially identical clustering across all three
# scales. In practice you can mitigate DBSCAN's scale sensitivity by
# standardizing your features (e.g. with
# :class:`~sklearn.preprocessing.StandardScaler`), but HDBSCAN removes
# the need to think about it at all.

# %%
# Summary
# -------
#
# - **DBSCAN** -- *Strength*: fast and simple when cluster density is
#   roughly uniform. *Limitation*: requires choosing ``eps``; breaks on
#   varying density and changes with data scale.
# - **OPTICS** -- *Strength*: the reachability plot exposes the full
#   hierarchy of density-based clusterings, and DBSCAN-like cuts can be
#   extracted from a single fit. *Limitation*: slower than DBSCAN; the
#   Xi extraction has its own hyperparameters.
# - **HDBSCAN** -- *Strength*: no ``eps`` to tune; adapts to clusters of
#   varying density and is scale-invariant. *Limitation*: higher memory
#   cost and slower than DBSCAN; tuning still benefits from
#   ``min_cluster_size`` and ``min_samples``.
#
# As a rule of thumb, start with DBSCAN when clusters look similarly
# dense, reach for HDBSCAN when densities vary or when the scale of the
# data is uncertain, and use OPTICS when you want to *inspect* the
# hierarchy of clusterings rather than commit to a single one.
#
# See the :ref:`clustering` user guide for more detail on tuning and
# evaluation of these algorithms.
