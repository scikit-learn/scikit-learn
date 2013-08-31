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

from collections import defaultdict
from operator import itemgetter

import numpy as np
import pylab as pl
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import EAC, KMeans, MSTCluster
from sklearn import metrics
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances


##############################################################################
# Generate sample data
n_samples = 150
noisy_circles = datasets.make_circles(n_samples=n_samples, factor=.5,
                                      noise=.05)
X, y_true = noisy_circles
main_metric = metrics.adjusted_mutual_info_score


##############################################################################
# Run K-Means many times with different k-values and record the distribution
# of adjusted mutual information scores.

print("Running k-means many times")
num_iterations = 50
km_ami_means = []
km_ami_std = []
k_values = list(range(2, 30))
for k in k_values:
    predictions = (KMeans(n_clusters=k).fit(X).labels_
                   for i in range(num_iterations))
    ami_scores = [main_metric(y_true, labels)
                  for labels in predictions]
    km_ami_means.append(np.mean(ami_scores))
    km_ami_std.append(np.std(ami_scores))


# Example k-means
km_labels = KMeans(n_clusters=2).fit(X).labels_


##############################################################################
# Run MSTCluster many times with different threshold values

print("Running MSTCluster many times")
num_iterations = 50
mst_scores = defaultdict(list)
D = pairwise_distances(X, metric="euclidean", n_jobs=1)
d_min = np.min(D)
d_max = np.max(D)
for threshold in np.arange(d_min, d_max, (d_max - d_min) / 100.):
    predictions = MSTCluster(threshold=threshold).fit(D).labels_
    ami_score = main_metric(y_true, predictions)
    n_clusters = len(set(predictions))
    mst_scores[n_clusters].append(ami_score)
for key in mst_scores:
    mst_scores[key] = np.mean(mst_scores[key])
xx = sorted(mst_scores.items(), key=itemgetter(0))
mst_values, mst_ami_means = zip(*xx)



##############################################################################
# Compute EAC
print("Running EAC Algorithm")
# Create a final clustering, allowing us to set the threshold value ourselves.
final_clusterer = MSTCluster(threshold=0.9)
model = EAC(final_clusterer=final_clusterer, random_state=42).fit(X)
eac_labels = model.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(eac_labels)) - (1 if -1 in eac_labels else 0)
print('Estimated number of clusters: %d' % n_clusters_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(y_true, eac_labels))
print("Completeness: %0.3f" % metrics.completeness_score(y_true, eac_labels))
print("V-measure: %0.3f" % metrics.v_measure_score(y_true, eac_labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(y_true, eac_labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(y_true, eac_labels))


##############################################################################
# Plot results

pl.title('EAC comparison with K-means')
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

# Subplot showing distribution of scores

# Plot two examples of k-means
ax = pl.subplot(3, 1, 1)
ax.scatter(X[:, 0], X[:,1], color=colors[km_labels].tolist(), s=10)


# Plot EAC labels
ax = pl.subplot(3, 1, 2)
ax.scatter(X[:, 0], X[:,1], color=colors[eac_labels].tolist(), s=10)

# Plot distribution of scores (from main_metric)
ax = pl.subplot(3, 1, 3)
# k-means
ax.plot(k_values, km_ami_means)
ax.errorbar(k_values, km_ami_means, yerr=km_ami_std, fmt='ro', label='k-means')

# MST
ax.plot(mst_values, mst_ami_means)
ax.errorbar(mst_values, mst_ami_means, fmt='g*', label='MST')
score = main_metric(y_true, eac_labels)
ax.scatter([n_clusters_,], [score,], label='EAC', s=40)
ax.legend()




pl.show()
