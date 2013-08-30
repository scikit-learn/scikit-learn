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
import pylab as pl
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import EAC, KMeans, MSTCluster
from sklearn import metrics
from sklearn import datasets
from sklearn.preprocessing import StandardScaler


##############################################################################
# Generate sample data
n_samples = 150
noisy_circles = datasets.make_circles(n_samples=n_samples, factor=.5,
                                      noise=.05)
X, y_true = noisy_circles
main_metric = metrics.adjusted_mutual_info_score


##############################################################################
# Run K-Means many times with random k-values and record the distribution
# of adjusted mutual information scores.

print("Plotting many runs of k-means")
num_iterations = 50
ami_means = []
ami_std = []
k_values = list(range(2, 30))
for k in k_values:
    predictions = (KMeans(n_clusters=k).fit(X).labels_
                   for i in range(num_iterations))
    ami_scores = [main_metric(y_true, labels)
                  for labels in predictions]
    ami_means.append(np.mean(ami_scores))
    ami_std.append(np.std(ami_scores))


# Example k-means
km_labels = (KMeans(n_clusters=i+2).fit(X).labels_
             for i in range(4))



##############################################################################
# Compute EAC
print("Running EAC Algorithm")
final_clusterer = MSTCluster(threshold=0.9)
model = EAC(final_clusterer=final_clusterer).fit(X)
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
for i, cur_km_labels in enumerate(km_labels):
    ax = pl.subplot(2, 4, i+1)
    ax.scatter(X[:, 0], X[:,1], color=colors[cur_km_labels].tolist(), s=10)


# Plot distribution of scores (from main_metric)
ax = pl.subplot(2, 4, 5)
ax.plot(k_values, ami_means)
ax.errorbar(k_values, ami_means, yerr=ami_std, fmt='ro', label='k-means')
score = main_metric(y_true, eac_labels)
print("Score: {:.4f}, n_clusters={}".format(score, n_clusters_))
ax.scatter([n_clusters_,], [score,], label='EAC')
ax.legend()


# Plot EAC labels
ax = pl.subplot(2, 4, 6)
ax.scatter(X[:, 0], X[:,1], color=colors[eac_labels].tolist(), s=10)



pl.show()
