# -*- coding: utf-8 -*-
"""
===================================
Demo of DBSCAN clustering algorithm
===================================

Reference:
Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
In: Proceedings of the 2nd International Conference on Knowledge Discovery
and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 2006

"""
print __doc__

import numpy as np
from scipy.spatial import distance
from scikits.learn.cluster import DBSCAN
from scikits.learn import metrics
from scikits.learn.datasets.samples_generator import make_blobs


##############################################################################
# Generate sample data
centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=300, centers=centers, cluster_std=0.4)

##############################################################################
# Compute similarities
D = distance.squareform(distance.pdist(X))
S = 1 - (D / np.max(D))

##############################################################################
# Compute DBSCAN
db = DBSCAN().fit(S, eps=0.95, min_points=5)
core_points = db.core_points_
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print 'Estimated number of clusters: %d' % n_clusters_
print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels)
print "Completeness: %0.3f" % metrics.completeness_score(labels_true, labels)
print "V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels)

##############################################################################
# Plot result
import pylab as pl
from itertools import cycle

pl.close('all')
pl.figure(1)
pl.clf()

# Black removed and is used for noise instead.
colors = cycle('bgrcmybgrcmybgrcmybgrcmy')
for k, col in zip(set(labels), colors):
    if k == -1:
        # Black used for noise.
        col = 'k'
        markersize = 6
    class_members = [index[0] for index in np.argwhere(labels == k)]
    cluster_core_points = [index for index in core_points
                           if labels[index] == k]
    for index in class_members:
        x = X[index]
        if index in core_points and k != -1:
            markersize = 14
        else:
            markersize = 6
        pl.plot(x[0], x[1], 'o', markerfacecolor=col,
                markeredgecolor='k', markersize=markersize)

pl.title('Estimated number of clusters: %d' % n_clusters_)
pl.show()
