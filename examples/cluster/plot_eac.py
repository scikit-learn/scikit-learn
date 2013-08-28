# -*- coding: utf-8 -*-
"""
================================
Demo of EAC clustering algorithm
================================

Uses many iterations of k-means with random values for k to create a
"co-association matrix", where C[i][j] is the frequency of times that instances
i and j are clustered together. The MSTCluster algorithm is run on this 
co-association matrix to form the final clusters.

"""
print(__doc__)

import numpy as np
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import EAC, MSTCluster
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler


##############################################################################
# Generate sample data
centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.3,
                            random_state=0)

X = StandardScaler().fit_transform(X)

##############################################################################
# Compute EAC
final_clusterer = MSTCluster(threshold=0.5)
model = EAC(final_clusterer=final_clusterer).fit(X)
labels = model.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters: %d' % n_clusters_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, labels))


##############################################################################
# Plot resulting dendrogram with colored labels representing the final clusters
import pylab as pl
ax = pl.subplot(111)
s = 40
ax.scatter(X[:, 0], X[:, 1], c=labels_true, s=s*4)
ax.scatter(X[:, 0], X[:, 1], c=labels, s=s)

pl.title('Estimated number of clusters: %d' % n_clusters_)
pl.show()
