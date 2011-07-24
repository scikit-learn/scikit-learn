# -*- coding: utf-8 -*-
"""
=============================================================
Demo of K-means clustering algorithm with Vector Quantisation
=============================================================

"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import datasets
from scikits.learn import metrics
from scikits.learn.cluster import KMeans

# Sample data taken from the Iris dataset.
iris = datasets.load_iris()
# These two features provide good basic separation of Iris.
# Features (0, 3) is better, but this graph looks better as an example of VQ.
X = iris.data[:, (1, 3)]
Y = iris.target
# Step size of the mesh.
h = .02

# Compute a k-means model on the data.
km = KMeans().fit(X, k=3)
centroids = km.cluster_centers_
labels = km.labels_

# Number of clusters in labels.
n_clusters_ = len(set(labels))

# Print some descriptive statistics about the quality of the clusters.
print 'Estimated number of clusters: %d' % n_clusters_
print "Homogeneity: %0.3f" % metrics.homogeneity_score(Y, labels)
print "Completeness: %0.3f" % metrics.completeness_score(Y, labels)
print "V-measure: %0.3f" % metrics.v_measure_score(Y, labels)


# Plot the decision boundary. For that, we will asign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = km.transform(np.c_[xx.ravel(), yy.ravel()], method='vq')

# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.set_cmap(pl.cm.Paired)
pl.pcolormesh(xx, yy, Z)

# Plot also the training points
pl.scatter(X[:, 0], X[:, 1], c=Y, marker='o', s=72)
# Plot the centroids a bit larger, only as an X
pl.scatter(centroids[:, 0], centroids[:, 1],
           marker='x', s=169, linewidths=3,
           color='w')
pl.title('K-means clustering algorithm with Vector Quantisation\n'
         'Centroids are marked with white cross')
pl.axis('tight')
pl.show()
