"""
===================================
Demo of OPTICS clustering algorithm
===================================

Finds core samples of high density and expands clusters from them.
This example uses data that is generated so that the clusters have
different densities.

The clustering is first used in its automatic settings, which is the
:class:`sklearn.cluster.OPTICS` algorithm, and then setting specific
thresholds on the reachability, which corresponds to DBSCAN.

We can see that the different clusters of OPTICS can be recovered with
different choices of thresholds in DBSCAN.

"""

# Authors: Shane Grigsby <refuge@rocktalus.com>
#          Amy X. Zhang <axz@mit.edu>
# License: BSD 3 clause


from sklearn.cluster import OPTICS
import matplotlib.gridspec as gridspec


import numpy as np

import matplotlib.pyplot as plt

# Generate sample data

np.random.seed(0)
n_points_per_cluster = 250

C1 = [-5, -2] + .8 * np.random.randn(n_points_per_cluster, 2)
C2 = [4, -1] + .1 * np.random.randn(n_points_per_cluster, 2)
C3 = [1, -2] + .2 * np.random.randn(n_points_per_cluster, 2)
C4 = [-2, 3] + .3 * np.random.randn(n_points_per_cluster, 2)
C5 = [3, -2] + 1.6 * np.random.randn(n_points_per_cluster, 2)
C6 = [5, 6] + 2 * np.random.randn(n_points_per_cluster, 2)
X = np.vstack((C1, C2, C3, C4, C5, C6))

clust = OPTICS(min_samples=9, rejection_ratio=0.5)

# Run the fit
clust.fit(X)

_, labels_025 = clust.extract_dbscan(0.25)
_, labels_075 = clust.extract_dbscan(0.75)

space = np.arange(len(X))
reachability = clust.reachability_[clust.ordering_]
labels = clust.labels_[clust.ordering_]

plt.figure(figsize=(10, 7))
G = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(G[0, :])
ax2 = plt.subplot(G[1, 0])
ax3 = plt.subplot(G[1, 1])
ax4 = plt.subplot(G[1, 2])

# Reachability plot
color = ['g.', 'r.', 'b.', 'y.', 'c.']
for k, c in zip(range(0, 5), color):
    Xk = space[labels == k]
    Rk = reachability[labels == k]
    ax1.plot(Xk, Rk, c, alpha=0.3)
ax1.plot(space[labels == -1], reachability[labels == -1], 'k.', alpha=0.3)
ax1.plot(space, np.full_like(space, 0.75, dtype=float), 'k-', alpha=0.5)
ax1.plot(space, np.full_like(space, 0.25, dtype=float), 'k-.', alpha=0.5)
ax1.set_ylabel('Reachability (epsilon distance)')
ax1.set_title('Reachability Plot')

# OPTICS
color = ['g.', 'r.', 'b.', 'y.', 'c.']
for k, c in zip(range(0, 5), color):
    Xk = X[clust.labels_ == k]
    ax2.plot(Xk[:, 0], Xk[:, 1], c, alpha=0.3)
ax2.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], 'k+', alpha=0.1)
ax2.set_title('Automatic Clustering\nOPTICS')

# DBSCAN at 0.25
color = ['g', 'greenyellow', 'olive', 'r', 'b', 'c']
for k, c in zip(range(0, 6), color):
    Xk = X[labels_025 == k]
    ax3.plot(Xk[:, 0], Xk[:, 1], c, alpha=0.3, marker='.')
ax3.plot(X[labels_025 == -1, 0], X[labels_025 == -1, 1], 'k+', alpha=0.1)
ax3.set_title('Clustering at 0.25 epsilon cut\nDBSCAN')

# DBSCAN at 0.75
color = ['g.', 'm.', 'y.', 'c.']
for k, c in zip(range(0, 4), color):
    Xk = X[labels_075 == k]
    ax4.plot(Xk[:, 0], Xk[:, 1], c, alpha=0.3)
ax4.plot(X[labels_075 == -1, 0], X[labels_075 == -1, 1], 'k+', alpha=0.1)
ax4.set_title('Clustering at 0.75 epsilon cut\nDBSCAN')

plt.tight_layout()
plt.show()
