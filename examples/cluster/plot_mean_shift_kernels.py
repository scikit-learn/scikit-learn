"""
============================================
Mean-Shift clustering with different kernels
============================================

This example shows the different outcomes with different kernels.
The flat kernel simply wants to have as many points as possible located
in it, while the rbf kernel wants as many points possible near its center.
This means that the flat kernel will be more likely to try to cover up two
or more distinct clusters because it doesn't care if the dense cluster
centers are located at its periphery.

With real world data the rbf kernel is probably preferred. The cluster center
is more likely to end up directly above an area with high density. It should
also be more robust on datasets with few points and noise.

The code borrows heavily from the other Mean-Shift example, and a
blog post by the original author:
http://sociograph.blogspot.nl/2011/11/accessible-introduction-to-mean-shift.html

Reference:

Dorin Comaniciu and Peter Meer, "Mean Shift: A robust approach toward
feature space analysis". IEEE Transactions on Pattern Analysis and
Machine Intelligence. 2002. pp. 603-619.

Yizong Cheng, "Mean Shift, Mode Seeking, and Clustering". IEEE
Transactions on Pattern Analysis and Machine Intelligence. 1995. pp. 790-799.

"""
print(__doc__)

import numpy as np
from sklearn.cluster import MeanShift
from sklearn.datasets.samples_generator import make_blobs

###############################################################################
# Generate sample data
centers = [[2, 0], [-2, 0]]
X, _ = make_blobs(n_samples=10000, centers=centers, cluster_std=(1.0, 1.0))

###############################################################################
# Compute clustering with MeanShift

# The difference between the kernels is best noticable when the bandwidth
# parameter is chosen somewhat too large. The bandwidth estimation function
# returns a value around 1.44
bandwidth = 3

ms_rbf = MeanShift(bandwidth=bandwidth, bin_seeding=True, kernel='rbf',
                   gamma=2)
ms_rbf.fit(X)
labels_rbf = ms_rbf.labels_
cluster_centers_rbf = ms_rbf.cluster_centers_

labels_unique_rbf = np.unique(labels_rbf)
n_clusters_rbf = len(labels_unique_rbf)

print("number of estimated clusters using the rbf kernel: %d" % n_clusters_rbf)


ms_flat = MeanShift(bandwidth=bandwidth, bin_seeding=True, kernel='flat')
ms_flat.fit(X)
labels_flat = ms_flat.labels_
cluster_centers_flat = ms_flat.cluster_centers_

labels_unique_flat = np.unique(labels_flat)
n_clusters_flat = len(labels_unique_flat)

print("number of estimated clusters using the flat kernel: %d" %
      n_clusters_flat)

###############################################################################
# Plot result
import pylab as pl
from itertools import cycle

pl.figure(1)
pl.clf()
colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

pl.subplot(2, 1, 1)
for k, col in zip(range(n_clusters_rbf), colors):
    my_members = labels_rbf == k
    cluster_center = cluster_centers_rbf[k]
    pl.plot(X[my_members, 0], X[my_members, 1], col + '.')
    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            markeredgecolor='k', markersize=14)
pl.title('Clustering using the rbf kernel')

pl.subplot(2, 1, 2)
for k, col in zip(range(n_clusters_flat), colors):
    my_members = labels_flat == k
    cluster_center = cluster_centers_flat[k]
    pl.plot(X[my_members, 0], X[my_members, 1], col + '.')
    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            markeredgecolor='k', markersize=14)
pl.title('Clustering using the flat kernel')

pl.show()
pl.close()
