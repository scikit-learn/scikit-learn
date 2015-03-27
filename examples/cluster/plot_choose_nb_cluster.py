"""
============================================
Selecting number of clusters on toy datasets
============================================

This example shows several algorithms to choose the number of clusters,
for a particular clustering algorithm on a particular dataset. It mainly
illustrates that some algorithms are faster, some algorithms only understand
convex clusters (first dataset) and some algorithms understand non-convex
clusters (second and third datasets).

The running times only give intuition of which algorithm is faster. Running time
highly depends on a datasets number of samples and number of features.
"""
print(__doc__)

import time
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster.spectral import SpectralClustering
from sklearn.datasets import make_blobs, make_moons, make_circles
from sklearn.metrics.cluster.calinski_harabaz_index import max_CH_index
from sklearn.metrics.cluster.stability import stability
from sklearn.metrics.cluster.distortion_jump import distortion_jump
from sklearn.metrics.cluster.gap_statistic import gap_statistic
from sklearn.metrics.cluster.unsupervised import silhouette_score
from sklearn.preprocessing import StandardScaler

n_samples = 1500
seed = 1
datasets = [
    make_blobs(n_samples=n_samples, random_state=seed),
    make_circles(n_samples=n_samples, factor=.5, noise=.05, shuffle=True, random_state=seed),
    make_moons(n_samples=n_samples, noise=.05, shuffle=True, random_state=seed),
]

cluster_estimator = SpectralClustering(eigen_solver='arpack', affinity="nearest_neighbors")


def max_silhouette(X, cluster_estimator, k_max=None):
    if not k_max:
        k_max = int(X.shape[0] / 2)
    silhouettes = []
    for k in range(2, k_max + 1):
        cluster_estimator.set_params(n_clusters=k)
        labels = cluster_estimator.fit_predict(X)
        silhouettes.append((k, silhouette_score(X, labels)))
    return max(silhouettes, key=itemgetter(1))[0]


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
nb_colors = len(colors)

plt.figure(figsize=(13, 9.5))
plt.subplots_adjust(left=.001, right=.999, bottom=.001, top=.96, wspace=.05,
                    hspace=.01)

plot_num = 1
printed_header = False

for dataset in datasets:
    X, true_labels = dataset
    # normalize dataset for nicer plotting
    X = StandardScaler().fit_transform(X)

    for name, func_choose_nb_cluster in {
            'Silhouette': max_silhouette,
            'Stability': stability,
            'Gap statistic': gap_statistic,
            'Calinski-Harabasz index': max_CH_index,
            'Distortion jump': distortion_jump,
    }.items():
        # predict cluster memberships
        t0 = time.time()
        nb_cluster = func_choose_nb_cluster(X, cluster_estimator, k_max=10)
        t1 = time.time()

        # retrieving clustering done
        cluster_estimator.set_params(n_clusters=nb_cluster)
        y_pred = cluster_estimator.fit_predict(X)

        # plot
        plt.subplot(3, 5, plot_num)
        if not printed_header:
            plt.title(name, size=18)
        points_color = [colors[y % nb_colors] for y in y_pred]
        plt.scatter(X[:, 0], X[:, 1], color=points_color, s=10)

        plt.xlim(-2, 2)
        plt.ylim(-2, 2)
        plt.xticks(())
        plt.yticks(())
        plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
                 transform=plt.gca().transAxes, size=15,
                 horizontalalignment='right')
        plot_num += 1
    printed_header = True

plt.show()
