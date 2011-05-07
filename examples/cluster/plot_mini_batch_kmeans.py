"""
=====================================================
A demo of the MiniBatch K Means clustering algorithm
=====================================================

"""
print __doc__

import numpy as np
from scikits.learn.cluster import MiniBatchKMeans

##############################################################################
# Generate sample data
np.random.seed(0)

n_points_per_cluster = 250
n_clusters = 3
n_points = n_points_per_cluster * n_clusters
means = np.array([[1, 1], [-1, -1], [1, -1]])
std = .6
iterations = 100
batch_size = 50

X = np.empty((0, 2))
for i in range(n_clusters):
    X = np.r_[X, means[i] + std * np.random.randn(n_points_per_cluster, 2)]

##############################################################################
# Compute clustering with MiniBatchKMeans
mbk = MiniBatchKMeans(init='k-means++',
                      k=3,
                      chunk_size=batch_size)
mbk.fit(X)
labels = mbk.labels_
cluster_centers = mbk.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

##############################################################################
# Plot result
import pylab as pl
from itertools import cycle

pl.figure(1)
pl.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    pl.plot(X[my_members, 0], X[my_members, 1], col + '.')
    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                                    markeredgecolor='k', markersize=14)
pl.title('Clustering with MiniBatchKMeans')
pl.show()
