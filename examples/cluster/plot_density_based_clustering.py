"""
===========================================
Demo of DBSCAN, OPTICS and HDBSCAN clustering
===========================================

This example compares three density-based clustering algorithms on the same
dataset: DBSCAN, OPTICS, and HDBSCAN. All three algorithms discover clusters
based on the density of points rather than their distance from a centroid,
which allows them to find clusters of arbitrary shapes and handle noise.

While the algorithms share the same core idea, they differ in how they handle
clusters of varying density. DBSCAN requires a fixed neighbourhood radius
(``eps``), which makes it sensitive to clusters of different densities. OPTICS
generalises DBSCAN by ordering points according to their reachability distance,
making the density structure of the data visible without committing to a single
``eps``. HDBSCAN extends this further by building a full cluster hierarchy and
extracting the most stable clusters automatically, removing the need to manually
inspect a reachability plot.

The dataset used here contains four clusters with intentionally different
densities to highlight where each algorithm succeeds and where it struggles.
"""
# %%
# Generating the dataset
# ----------------------
#
# We use :func:`~sklearn.datasets.make_blobs` to generate a dataset with four
# clusters of intentionally different densities. Two clusters are tight
# (small ``cluster_std``) and two are spread out (large ``cluster_std``).
# This mix of densities is what makes the comparison between the three
# algorithms meaningful.

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
# Visualising the ground truth
# ----------------------------
#
# Before running any algorithm, let's look at the raw dataset and its true
# cluster labels. This gives us a baseline to compare each algorithm against.

fig, ax = plt.subplots(figsize=(8, 5))
unique_labels = set(labels_true)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

for k, col in zip(unique_labels, colors):
    mask = labels_true == k
    ax.scatter(
        X[mask, 0],
        X[mask, 1],
        c=[col],
        s=14,
        label=f"Cluster {k}",
    )

ax.set_title("Ground truth clusters")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# DBSCAN: simple but sensitive to density differences
# ----------------------------------------------------
#
# DBSCAN groups together points that are closely packed by requiring that a
# minimum number of points (``min_samples``) lie within a fixed radius
# (``eps``) of each other. Points that lie in low-density regions are labelled
# as noise (-1).
#
# The key limitation of DBSCAN is that ``eps`` is a global threshold applied
# uniformly to the entire dataset. When clusters have different densities —
# as in our dataset — no single value of ``eps`` works perfectly for all of
# them. A value tuned for the tight clusters will fragment the sparse ones
# into noise, while a value tuned for the sparse clusters will merge the
# tight ones together.
#
# We demonstrate this below by running DBSCAN with two different ``eps``
# values on the same dataset.


def plot_dbscan(X, labels, eps, ax):
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    n_noise = list(labels).count(-1)

    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]
        mask = labels == k
        marker = "o" if k != -1 else "x"
        ax.scatter(X[mask, 0], X[mask, 1], c=[col], s=14, marker=marker)

    ax.set_title(
        f"DBSCAN: eps={eps}\n"
        f"{n_clusters} clusters, {n_noise} noise points"
    )


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# eps tuned for tight clusters — sparse clusters fragment into noise
db1 = DBSCAN(eps=0.15, min_samples=5).fit(X)
plot_dbscan(X, db1.labels_, eps=0.15, ax=ax1)

# eps tuned for sparse clusters — tight clusters merge incorrectly
db2 = DBSCAN(eps=0.6, min_samples=5).fit(X)
plot_dbscan(X, db2.labels_, eps=0.6, ax=ax2)

fig.suptitle(
    "DBSCAN: no single eps works for clusters of different density",
    fontsize=13,
)
plt.tight_layout()
plt.show()

# %%
# The plots above illustrate the core weakness of DBSCAN. With ``eps=0.15``
# (left), the algorithm correctly separates the two tight clusters but treats
# the two sparse clusters entirely as noise. With ``eps=0.6`` (right), the
# sparse clusters are recovered, but the tight clusters are incorrectly merged
# into a single large cluster.
#
# This is not a failure of DBSCAN — it is a fundamental limitation of using
# a single global density threshold. OPTICS was designed specifically to
# address this problem.

# %%
# OPTICS: revealing density structure without committing to epsilon
# ----------------------------------------------------------------
#
# OPTICS (Ordering Points To Identify the Clustering Structure) generalises
# DBSCAN by computing a reachability distance for every point — a measure of
# how hard it is to "reach" a point from a dense cluster core. It then orders
# all points so that nearby points in dense regions appear close together in
# the ordering.
#
# The result is a reachability plot: a bar chart where valleys correspond to
# dense clusters and peaks correspond to sparse regions or noise. Crucially,
# clusters of different densities appear as valleys at different heights,
# making the full density structure of the data visible in a single plot.
#
# Once OPTICS has been fitted, you can extract flat cluster labels at any
# epsilon threshold using :func:`~sklearn.cluster.cluster_optics_dbscan`,
# effectively running DBSCAN at multiple thresholds without re-fitting.

clust = OPTICS(min_samples=10, xi=0.05, min_cluster_size=0.05)
clust.fit(X)

# Extract labels at two different epsilon thresholds
labels_200 = cluster_optics_dbscan(
    reachability=clust.reachability_,
    core_distances=clust.core_distances_,
    ordering=clust.ordering_,
    eps=0.2,
)
labels_075 = cluster_optics_dbscan(
    reachability=clust.reachability_,
    core_distances=clust.core_distances_,
    ordering=clust.ordering_,
    eps=0.75,
)

# %%
# Reachability plot
# ~~~~~~~~~~~~~~~~~
#
# The reachability plot below shows the density structure OPTICS discovered.
# Each bar represents one point, ordered by the OPTICS traversal. Valleys
# (short bars) are dense regions — clusters. Peaks (tall bars) are sparse
# regions or transitions between clusters.
#
# The dashed horizontal lines show two epsilon thresholds. Reading off the
# clusters at each threshold gives us the flat clusterings shown afterwards.

reachability = clust.reachability_[clust.ordering_]
labels_optics = clust.labels_[clust.ordering_]

fig, ax = plt.subplots(figsize=(10, 4))
colors = ["g", "steelblue", "orange", "purple"]

for klass, color in zip(range(0, 4), colors):
    Xk = reachability[labels_optics == klass]
    ax.plot(
        np.where(labels_optics == klass)[0],
        Xk,
        color,
        alpha=0.7,
        linewidth=0.8,
    )

ax.plot(
    np.where(labels_optics == -1)[0],
    reachability[labels_optics == -1],
    "k.",
    alpha=0.3,
    markersize=2,
)
ax.axhline(0.2, linestyle="--", color="k", linewidth=1, label="eps=0.20")
ax.axhline(0.75, linestyle="--", color="r", linewidth=1, label="eps=0.75")
ax.set_ylabel("Reachability distance")
ax.set_xlabel("Point ordering")
ax.set_title("OPTICS reachability plot\nValleys = clusters, peaks = noise or boundaries")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# Cluster assignments from OPTICS at two epsilon thresholds
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The reachability plot gives us the full picture. We can now extract flat
# cluster labels at any epsilon threshold without re-running the algorithm.
# Below we compare the cluster assignments at ``eps=0.20`` and ``eps=0.75``.


def plot_optics_labels(X, labels, title, ax):
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    n_noise = list(labels).count(-1)

    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]
        mask = labels == k
        marker = "o" if k != -1 else "x"
        ax.scatter(X[mask, 0], X[mask, 1], c=[col], s=14, marker=marker)

    ax.set_title(f"{title}\n{n_clusters} clusters, {n_noise} noise points")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plot_optics_labels(X, labels_200, "OPTICS at eps=0.20", ax1)
plot_optics_labels(X, labels_075, "OPTICS at eps=0.75", ax2)
fig.suptitle(
    "OPTICS: same fit, different epsilon thresholds",
    fontsize=13,
)
plt.tight_layout()
plt.show()

# %%
# OPTICS successfully captures the density structure of the data, but
# extracting good clusters still requires the user to inspect the reachability
# plot and choose an appropriate threshold. HDBSCAN automates this final step.

# %%
# HDBSCAN: automatic multi-scale clustering
# -----------------------------------------
#
# HDBSCAN (Hierarchical DBSCAN) builds on the same mutual reachability
# distances as OPTICS but goes one step further: it constructs a full
# cluster hierarchy across all density levels and then automatically extracts
# the clusters that are most stable — those that persist across the widest
# range of density thresholds before splitting or dissolving into noise.
#
# This means HDBSCAN requires no ``eps`` parameter at all. The only required
# parameter is ``min_cluster_size``, which has an intuitive meaning: the
# smallest grouping of points you consider a meaningful cluster.
#
# A unique feature of HDBSCAN is its ``probabilities_`` attribute. Every
# point receives a membership probability between 0 and 1 indicating how
# confidently it belongs to its assigned cluster. Points deep inside a dense
# core score close to 1.0, while points on the fuzzy boundary of a cluster
# score lower. In the plot below, marker size reflects this probability —
# smaller markers are less certain assignments.

hdb = HDBSCAN(min_cluster_size=15, min_samples=5, copy=True)
hdb.fit(X)

fig, ax = plt.subplots(figsize=(8, 5))
unique_labels = set(hdb.labels_)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
n_noise = list(hdb.labels_).count(-1)

for k, col in zip(unique_labels, colors):
    if k == -1:
        col = [0, 0, 0, 1]
    mask = hdb.labels_ == k
    # marker size reflects membership probability
    sizes = 5 + 40 * hdb.probabilities_[mask]
    marker = "o" if k != -1 else "x"
    ax.scatter(
        X[mask, 0],
        X[mask, 1],
        c=[col],
        s=sizes if k != -1 else 14,
        marker=marker,
        alpha=0.8,
    )

ax.set_title(
    f"HDBSCAN: min_cluster_size=15\n"
    f"{n_clusters} clusters, {n_noise} noise points\n"
    f"Marker size = membership probability"
)
plt.tight_layout()
plt.show()

# %%
# Side-by-side comparison
# -----------------------
#
# Finally, we compare all three algorithms on the same dataset in a single
# figure. DBSCAN is shown with the best manually tuned ``eps`` value.
# OPTICS is shown at the threshold that best recovers all four clusters.
# HDBSCAN requires no threshold tuning at all.
#
# This comparison highlights the key trade-offs:
#
# - **DBSCAN** is the simplest and fastest, but requires careful tuning of
#   ``eps`` and struggles when clusters have different densities.
# - **OPTICS** reveals the full density structure of the data through the
#   reachability plot, making it easier to choose a good threshold, but still
#   requires human judgment to extract final clusters.
# - **HDBSCAN** handles varying densities automatically and provides
#   membership probabilities, making it the most robust choice for real-world
#   data where cluster densities are unknown in advance.

db_best = DBSCAN(eps=0.3, min_samples=5).fit(X)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
titles = ["DBSCAN (eps=0.30)", "OPTICS (eps=0.75)", "HDBSCAN (min_cluster_size=15)"]
all_labels = [db_best.labels_, labels_075, hdb.labels_]

for ax, labels, title in zip(axes, all_labels, titles):
    unique_labels = set(labels)
    colors = [
        plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))
    ]
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    n_noise = list(labels).count(-1)

    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]
        mask = labels == k
        marker = "o" if k != -1 else "x"
        ax.scatter(
            X[mask, 0],
            X[mask, 1],
            c=[col],
            s=14,
            marker=marker,
            alpha=0.8,
        )
    ax.set_title(f"{title}\n{n_clusters} clusters, {n_noise} noise points")

fig.suptitle(
    "Comparing DBSCAN, OPTICS, and HDBSCAN on a dataset with varying cluster densities",
    fontsize=12,
)
plt.tight_layout()
plt.show()