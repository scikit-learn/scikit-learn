"""
===================================
Demo of OPTICS clustering algorithm
===================================

Finds core samples of high density and expands clusters from them.
The key specificity of this example is that the data is generated
so that the clusters have different densities.
"""

# Authors: Shane Grigsby <refuge@rocktalus.com>
#          Amy X. Zhang <axz@mit.edu>
# License: BSD 3 clause


from sklearn.cluster.optics_ import OPTICS
import numpy as np

import matplotlib.pyplot as plt

# Generate sample data

np.random.seed(0)
n_points_per_cluster = 250

X = np.empty((0, 2))
X = np.r_[X, [-5, -2] + .8 * np.random.randn(n_points_per_cluster, 2)]
X = np.r_[X, [4, -1] + .1 * np.random.randn(n_points_per_cluster, 2)]
X = np.r_[X, [1, -2] + .2 * np.random.randn(n_points_per_cluster, 2)]
X = np.r_[X, [-2, 3] + .3 * np.random.randn(n_points_per_cluster, 2)]
X = np.r_[X, [3, -2] + 1.6 * np.random.randn(n_points_per_cluster, 2)]
X = np.r_[X, [5, 6] + 2 * np.random.randn(n_points_per_cluster, 2)]

# Plot scatterplot of points

plt.figure(figsize=(7, 7))

plt.subplot(221)
plt.plot(X[:, 0], X[:, 1], 'b.', ms=2)
plt.title("Raw Data")

# Compute OPTICS

clust = OPTICS(eps=30.3, min_samples=9, metric='minkowski')

# Run the fit
clust.fit(X)

# Plot result

core_samples_mask = np.zeros_like(clust.labels_, dtype=bool)
core_samples_mask[clust.core_sample_indices_] = True

# Black removed and is used for noise instead.
unique_labels = set(clust.labels_)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

plt.subplot(222)

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (clust.labels_ == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14, alpha=0.5)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=2, alpha=0.5)

plt.title("Automatic Clustering \n Estimated number of clusters: %d"
          % clust.n_clusters)

# (Re)-extract clustering structure, using a single eps to show comparison
# with DBSCAN. This can be run for any clustering distance, and can be run
# multiple times without rerunning OPTICS. OPTICS does need to be re-run to c
# hange the min-pts parameter.

clust.extract(.15, 'dbscan')

core_samples_mask = np.zeros_like(clust.labels_, dtype=bool)
core_samples_mask[clust.core_sample_indices_] = True

# Black removed and is used for noise instead.
unique_labels = set(clust.labels_)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

plt.subplot(223)

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (clust.labels_ == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], '.', markerfacecolor=col,
             markeredgecolor='k', markersize=14, alpha=0.5)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], '.', markerfacecolor=col,
             markeredgecolor='k', markersize=2, alpha=0.5)

plt.title('DBSCAN with epsilon of 0.15 \n Estimated number of clusters: %d'
          % clust.n_clusters)

# Try with different eps to highlight the problem


clust.extract(.4, 'dbscan')


core_samples_mask = np.zeros_like(clust.labels_, dtype=bool)
core_samples_mask[clust.core_sample_indices_] = True

# Black removed and is used for noise instead.
unique_labels = set(clust.labels_)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

plt.subplot(224)

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (clust.labels_ == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], '.', markerfacecolor=col,
             markeredgecolor='k', markersize=14, alpha=0.5)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], '.', markerfacecolor=col,
             markeredgecolor='k', markersize=2, alpha=0.5)

plt.title('DBSCAN with epsilon of 0.40 \n Estimated number of clusters: %d'
          % clust.n_clusters)

plt.tight_layout()
plt.show()
