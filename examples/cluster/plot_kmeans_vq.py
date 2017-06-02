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
from scikits.learn.cluster import KMeans
from scikits.learn import metrics
from scikits.learn.datasets.samples_generator import make_blobs


##############################################################################
# Generate sample data
centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=300, centers=centers, cluster_std=0.5)


##############################################################################
# Compute DBSCAN
km = KMeans().fit(X, k=3)
centroids = km.cluster_centers_
labels = km.labels_

# Number of clusters in labels
n_clusters_ = len(set(labels))

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

colors = 'bgrcmybgrcmybgrcmybgrcmy'

# Plot the background image
# make these smaller to increase the resolution
dx, dy = 0.1, 0.1
x_values = np.arange(-3.0, 3.0, dx)
y_values = np.arange(-3.0, 3.0, dy)
markersize = 9
for x in x_values:
    for y in y_values:
        point = np.array([[x, y],])
        label = km.transform(point, method='vq')[0]
        pl.plot(x, y, 'o', markerfacecolor=colors[label],
                markeredgecolor=colors[label], markersize=markersize)
        


# Overlay the original points
markersize = 6
for k, col in zip(set(labels), colors):
    class_members = [index[0] for index in np.argwhere(labels == k)]
    for index in class_members:
        x = X[index]
        pl.plot(x[0], x[1], 'o', markerfacecolor=col,
                markeredgecolor='k', markersize=markersize)

# Overlay the centroids
markersize = 12
for center in km.cluster_centers_:
    pl.plot(center[0], center[1], 'o', markerfacecolor='k',
            markeredgecolor='k', markersize=markersize)
pl.title('Estimated number of clusters: %d' % n_clusters_)
pl.show()
