"""
=====================================================
A demo of the K Means clustering algorithm
=====================================================

"""
print __doc__

import numpy as np
from scikits.learn.cluster import MiniBatchKMeans, KMeans

##############################################################################
# Generate sample data
np.random.seed(0)

n_points_per_cluster = 250
n_clusters = 3
n_points = n_points_per_cluster * n_clusters
means = np.array([[1, 1], [-1, -1], [1, -1]])
std = .6
batch_size = 50

X = np.empty((0, 2))
for i in range(n_clusters):
    X = np.r_[X, means[i] + std * np.random.randn(n_points_per_cluster, 2)]

##############################################################################
# Compute clustering with Means
k_means = KMeans(init='k-means++',
                 k=3)
k_means.fit(X)
k_means_labels = k_means.labels_
k_means_cluster_centers = k_means.cluster_centers_
k_means_labels_unique = np.unique(k_means_labels)


##############################################################################
# Compute clustering with MiniBatchKMeans
mbk = MiniBatchKMeans(init='k-means++',
                      k=3,
                      chunk_size=batch_size)
mbk.fit(X)
mbk_means_labels = mbk.labels_
mbk_means_cluster_centers = mbk.cluster_centers_
mbk_means_labels_unique = np.unique(mbk_means_labels)

##############################################################################
# Plot result
import matplotlib.pyplot as plt
from itertools import cycle

fig = plt.figure()
colors = cycle('bgr')

# KMeans
ax = fig.add_subplot(1, 3, 1)
for k, col in zip(range(n_clusters), colors):
    my_members = k_means_labels == k
    cluster_center = k_means_cluster_centers[k]
    ax.plot(X[my_members, 0], X[my_members, 1], col + '.')
    ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                                    markeredgecolor='k', markersize=14)
ax.set_title('Clustering with BatchKMeans')

# MiniBatchKMeans
ax = fig.add_subplot(1, 3, 2)
for k, col in zip(range(n_clusters), colors):
    my_members = mbk_means_labels == k
    cluster_center = mbk_means_cluster_centers[k]
    ax.plot(X[my_members, 0], X[my_members, 1], col + '.')
    ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                                    markeredgecolor='k', markersize=14)
ax.set_title('Clustering with MiniBatchKMeans')


identic = (mbk_means_labels == k_means_labels)
different = (mbk_means_labels != k_means_labels)


ax = fig.add_subplot(1, 3, 3)
for k, col in zip([identic, different], cycle('km')):
    my_members = k
    ax.plot(X[identic, 0], X[identic, 1], col + '.')

ax.set_title('')

plt.show()
