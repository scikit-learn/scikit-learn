# -*- coding: utf-8 -*-
"""
================================
Demo of EAC clustering algorithm
================================

Uses many iterations of k-means with random values for k to create a
"co-association matrix", where C[i][j] is the frequency of times that instances
i and j are clustered together. The SingleLinkageCluster algorithm is run on
this co-association matrix to form the final clusters.

"""
print(__doc__)

from collections import defaultdict
from operator import itemgetter

import numpy as np
import pylab as pl
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import eac, EAC, KMeans, SingleLinkageCluster
from sklearn import metrics
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances


##############################################################################
# Generate sample data
n_samples = 750
noisy_circles = datasets.make_circles(n_samples=n_samples, factor=.5,
                                      noise=.05)
X, y_true = noisy_circles
main_metric = metrics.adjusted_mutual_info_score



##############################################################################
# Compute EAC
print("Running EAC Algorithm")
# Create a final clustering, allowing us to set the threshold value ourselves.
threshold = 0.8
final_clusterer = SingleLinkageCluster(threshold=threshold,
                                       metric='precomputed')
model = EAC(default_final_clusterer=final_clusterer, random_state=42).fit(X)
y_pred = model.labels_
span_tree = model.final_clusterer.span_tree

print("Threshold: {:.4f}".format(threshold))
# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(y_pred)) - (1 if -1 in y_pred else 0)
print('Estimated number of clusters: %d' % n_clusters_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(y_true, y_pred))
print("Completeness: %0.3f" % metrics.completeness_score(y_true, y_pred))
print("V-measure: %0.3f" % metrics.v_measure_score(y_true, y_pred))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(y_true, y_pred))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(y_true, y_pred))
print("")

##############################################################################
# Run K-Means many times with different k-values and record the distribution
# of adjusted mutual information scores.

print("Running k-means many times")
num_iterations = 10
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
# Run SingleLinkageCluster many times with different threshold values

print("Running SingleLinkageCluster many times")
num_iterations = 50
sl_scores = defaultdict(list)
D = pairwise_distances(X, metric="euclidean", n_jobs=1)
d_min = np.min(D)
d_max = np.max(D)
for threshold in np.arange(d_min, d_max, (d_max - d_min) / 100.):
    model = SingleLinkageCluster(threshold=threshold, metric='euclidean')
    predictions = model.fit(D).labels_
    score = main_metric(y_true, predictions)
    n_clusters = len(set(predictions))
    sl_scores[n_clusters].append(score)
for key in sl_scores:
    sl_scores[key] = np.mean(sl_scores[key])
xx = sorted(sl_scores.items(), key=itemgetter(0))
sl_values, sl_ami_means = zip(*xx)


##############################################################################
# Plot results

fig = pl.figure(figsize=(8, 12))
fig.suptitle('EAC comparison with K-means')

colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 100)
# Subplot showing distribution of scores

# Plot two examples of k-means
ax = fig.add_subplot(3, 1, 1)
ax.scatter(X[:, 0], X[:,1], color=colors[km_labels].tolist(), s=10)
ax.set_title("K-means")

# Plot EAC labels
ax = fig.add_subplot(3, 1, 2)
ax.scatter(X[:, 0], X[:,1], color=colors[y_pred].tolist(), s=10)
ax.set_title("EAC")

# Plot distribution of scores (from main_metric)
ax = fig.add_subplot(3, 1, 3)
# k-means
ax.plot(k_values, km_ami_means)
ax.errorbar(k_values, km_ami_means, yerr=km_ami_std, fmt='ro', label='k-means')

# SingleLinkageCluster
ax.plot(sl_values, sl_ami_means)
ax.errorbar(sl_values, sl_ami_means, fmt='g*', label='Single Linkage')
score = main_metric(y_true, y_pred)
ax.scatter([n_clusters_,], [score,], label='EAC', s=40)
ax.legend()
ax.set_title("V-measure comparison")



pl.show()
