"""
=================================================
Demo of affinity propagation clustering algorithm
=================================================

Reference:
Brendan J. Frey and Delbert Dueck, "Clustering by Passing Messages
Between Data Points", Science Feb. 2007

"""
print __doc__

import numpy as np
from scikits.learn.cluster import AffinityPropagation
from scikits.learn.metrics import v_measure_score
from scikits.learn.metrics import homogeneity_score
from scikits.learn.metrics import completeness_score

##############################################################################
# Generate sample data
##############################################################################
np.random.seed(0)

means = (
    [1, 1],
    [-1, -1],
    [1, -1]
)
std = .5
n_clusters = len(means)
n_samples_per_cluster = 100

clusters = []
labels_true = []
for i, mean in enumerate(means):
    clusters.append(np.random.normal(
        loc=mean, scale=std, size=(n_samples_per_cluster, 2)))
    labels_true += [i] * n_samples_per_cluster

X = np.concatenate(clusters)

##############################################################################
# Compute similarities
##############################################################################
X_norms = np.sum(X * X, axis=1)
S = - X_norms[:, np.newaxis] - X_norms[np.newaxis, :] + 2 * np.dot(X, X.T)
p = 10 * np.median(S)

##############################################################################
# Compute Affinity Propagation
##############################################################################

af = AffinityPropagation().fit(S, p)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

n_clusters_ = len(cluster_centers_indices)

print 'Estimated number of clusters: %d' % n_clusters_
print "Homogeneity: %0.3f" % homogeneity_score(labels_true, labels)
print "Completeness: %0.3f" % completeness_score(labels_true, labels)
print "V-measure: %0.3f" % v_measure_score(labels_true, labels)


##############################################################################
# Plot result
##############################################################################

import pylab as pl
from itertools import cycle

pl.close('all')
pl.figure(1)
pl.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    class_members = labels == k
    cluster_center = X[cluster_centers_indices[k]]
    pl.plot(X[class_members, 0], X[class_members, 1], col + '.')
    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            markeredgecolor='k', markersize=14)
    for x in X[class_members]:
        pl.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

pl.title('Estimated number of clusters: %d' % n_clusters_)
pl.show()
