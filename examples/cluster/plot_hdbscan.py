# -*- coding: utf-8 -*-
"""
===================================
Demo of HDBSCAN clustering algorithm
===================================

Finds a clustering that has the greatest stability over a range
of epsilon values for standard DBSCAN. This allows clusterings
of different densities unlike DBSCAN.

"""
print(__doc__)

import numpy as np

from sklearn.cluster import HDBSCAN, DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import time

def make_var_density_blobs(n_samples=750, centers=[[0,0]], cluster_std=[0.5], random_state=0):
    samples_per_blob = n_samples // len(centers)
    blobs = [make_blobs(n_samples=samples_per_blob, centers=[c], cluster_std=cluster_std[i])[0]
             for i, c in enumerate(centers)]
    labels = [i * np.ones(samples_per_blob) for i in range(len(centers))]
    return np.vstack(blobs), np.hstack(labels)
        

##############################################################################
# Generate sample data
centers = [[1, 1], [-1, -1], [1, -1]]
densities = [0.2, 0.35, 0.5]
X, labels_true = make_var_density_blobs(n_samples=750, centers=centers, cluster_std=densities,
                            random_state=0)

X = StandardScaler().fit_transform(X)

##############################################################################
# Compute DBSCAN
hdb_t1 = time.time()
hdb = HDBSCAN(min_cluster_size=10).fit(X)
hdb_labels = hdb.labels_
hdb_elapsed_time = time.time() - hdb_t1

db_t1 = time.time()
db = DBSCAN().fit(X)
db_labels = db.labels_
db_elapsed_time = time.time() - db_t1

# Number of clusters in labels, ignoring noise if present.
n_clusters_hdb_ = len(set(hdb_labels)) - (1 if -1 in hdb_labels else 0)

print('\n\n++ HDBSCAN Results')
print('Estimated number of clusters: %d' % n_clusters_hdb_)
print('Elapsed time to cluster: %.4f s' % hdb_elapsed_time)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, hdb_labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, hdb_labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, hdb_labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, hdb_labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, hdb_labels))
print("Silhouette Coefficient: %0.3f"
      % metrics.silhouette_score(X, hdb_labels))

n_clusters_db_ = len(set(db_labels)) - (1 if -1 in db_labels else 0)

print('\n\n++ DBSCAN Results')
print('Estimated number of clusters: %d' % n_clusters_db_)
print('Elapsed time to cluster: %.4f s' % db_elapsed_time)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, db_labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, db_labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, db_labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, db_labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, db_labels))
print("Silhouette Coefficient: %0.3f"
      % metrics.silhouette_score(X, db_labels))

##############################################################################
# Plot result
import matplotlib.pyplot as plt

# Black removed and is used for noise instead.
hdb_unique_labels = set(hdb_labels)
db_unique_labels = set(db_labels)
hdb_colors = plt.cm.Spectral(np.linspace(0, 1, len(hdb_unique_labels)))
db_colors = plt.cm.Spectral(np.linspace(0, 1, len(db_unique_labels)))
fig = plt.figure()
hdb_axis = fig.add_subplot('121')
db_axis = fig.add_subplot('122')
for k, col in zip(hdb_unique_labels, hdb_colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    hdb_axis.plot(X[hdb_labels == k, 0], X[hdb_labels == k, 1], 'o', markerfacecolor=col,
                  markeredgecolor='k', markersize=6)
for k, col in zip(db_unique_labels, db_colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    db_axis.plot(X[db_labels == k, 0], X[db_labels == k, 1], 'o', markerfacecolor=col,
                  markeredgecolor='k', markersize=6)

hdb_axis.set_title('HDBSCAN\nEstimated number of clusters: %d' % n_clusters_hdb_)
db_axis.set_title('DBSCAN\nEstimated number of clusters: %d' % n_clusters_db_)
plt.show()
