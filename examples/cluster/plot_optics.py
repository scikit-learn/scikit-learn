"""
===================================
Demo of OPTICS clustering algorithm
===================================

Finds core samples of high density and expands clusters from them.
This example uses data that is generated so that the clusters have
different densities.
"""

# Authors: Shane Grigsby <refuge@rocktalus.com>
#          Amy X. Zhang <axz@mit.edu>
# License: BSD 3 clause


from sklearn.cluster.optics_ import OPTICS
import matplotlib.gridspec as gridspec


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


clust = OPTICS(min_samples=9, rejection_ratio=.5)

# Run the fit
clust.fit(X)

_, labels_025 = clust.extract_dbscan(0.25)
_, labels_075 = clust.extract_dbscan(0.75)

space = np.r_[0:len(X):1]
R = clust.reachability_[clust.ordering_]
C = clust.labels_[clust.ordering_]

plt.figure(figsize=(10, 7))
G = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(G[0, :])
ax1.set_ylabel('Reachability (epsilon distance)')
ax1.set_title('Reachability Plot')
ax2 = plt.subplot(G[1, 0])
ax2.set_title('Automatic Clustering')
ax3 = plt.subplot(G[1, 1])
ax3.set_title('Clustering at 0.25 epsilon cut')
ax4 = plt.subplot(G[1, 2])
ax4.set_title('Clustering at 0.75 epsilon cut')

# Reachability plot
ax1.plot(space[C == 1], R[C == 1], 'g.', alpha=0.3)
ax1.plot(space[C == 2], R[C == 2], 'r.', alpha=0.3)
ax1.plot(space[C == 3], R[C == 3], 'b.', alpha=0.3)
ax1.plot(space[C == 4], R[C == 4], 'y.', alpha=0.3)
ax1.plot(space[C == 5], R[C == 5], 'c.', alpha=0.3)
ax1.plot(space[C == -1], R[C == -1], 'k.', alpha=0.3)
ax1.plot(space, np.ones_like(space) * 0.75, 'k-', alpha=0.5)
ax1.plot(space, np.ones_like(space) * 0.25, 'k-.', alpha=0.5)

# OPTICS
ax2.plot(X[clust.labels_ == 1, 0], X[clust.labels_ == 1, 1], 'g.', alpha=0.3)
ax2.plot(X[clust.labels_ == 2, 0], X[clust.labels_ == 2, 1], 'r.', alpha=0.3)
ax2.plot(X[clust.labels_ == 3, 0], X[clust.labels_ == 3, 1], 'b.', alpha=0.3)
ax2.plot(X[clust.labels_ == 4, 0], X[clust.labels_ == 4, 1], 'y.', alpha=0.3)
ax2.plot(X[clust.labels_ == 5, 0], X[clust.labels_ == 5, 1], 'c.', alpha=0.3)
ax2.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], 'k+', alpha=0.1)

# DBSCAN at 0.25
ax3.plot(X[labels_025 == 1, 0], X[labels_025 == 1, 1], 'g.', alpha=0.3)
ax3.plot(X[labels_025 == 2, 0], X[labels_025 == 2, 1], 'm.', alpha=0.3)
ax3.plot(X[labels_025 == 3, 0], X[labels_025 == 3, 1], 'k.', alpha=0.3)
ax3.plot(X[labels_025 == 4, 0], X[labels_025 == 4, 1], 'r.', alpha=0.3)
ax3.plot(X[labels_025 == 5, 0], X[labels_025 == 5, 1], 'b.', alpha=0.3)
ax3.plot(X[labels_025 == 6, 0], X[labels_025 == 6, 1], 'c.', alpha=0.3)
ax3.plot(X[labels_025 == -1, 0], X[labels_025 == -1, 1], 'k+', alpha=0.1)

ax4.plot(X[labels_075 == 1, 0], X[labels_075 == 1, 1], 'g.', alpha=0.3)
ax4.plot(X[labels_075 == 2, 0], X[labels_075 == 2, 1], 'm.', alpha=0.3)
ax4.plot(X[labels_075 == 3, 0], X[labels_075 == 3, 1], 'y.', alpha=0.3)
ax4.plot(X[labels_075 == 4, 0], X[labels_075 == 4, 1], 'c.', alpha=0.3)
ax4.plot(X[labels_075 == -1, 0], X[labels_075 == -1, 1], 'k+', alpha=0.1)

plt.tight_layout()
plt.show()
