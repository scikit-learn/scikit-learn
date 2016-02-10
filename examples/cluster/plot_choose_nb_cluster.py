from __future__ import division
import time
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np

from sklearn.preprocessing import StandardScaler

# ----- importing datasets ------
from sklearn.datasets import make_blobs, make_moons, make_circles
n_samples = 1500
seed = 1
datasets = [
    make_blobs(n_samples=n_samples, random_state=seed),
    make_circles(n_samples=n_samples, factor=.5, noise=.05, shuffle=True, random_state=seed),
    make_moons(n_samples=n_samples, noise=.05, shuffle=True, random_state=seed),
]

# ----- Cluster estimator -----
from sklearn.cluster.spectral import SpectralClustering
cluster_estimator = SpectralClustering(eigen_solver='arpack', affinity="nearest_neighbors")


# ----- nb cluster estimators ------
from sklearn.metrics.cluster.calinski_harabaz_index import max_CH_index
from sklearn.metrics.cluster.stability import stability
from sklearn.metrics.cluster.distortion_jump import distortion_jump
from sklearn.metrics.cluster.gap_statistic import gap_statistic
from sklearn.metrics.cluster.unsupervised import silhouette_score


def max_silhouette(X, cluster_estimator, k_max=None):
    if not k_max:
        k_max = X.shape[0] // 2
    silhouettes = []
    for k in range(2, k_max + 1):
        cluster_estimator.set_params(n_clusters=k)
        labels = cluster_estimator.fit_predict(X)
        silhouettes.append((k, silhouette_score(X, labels)))
    return max(silhouettes, key=itemgetter(1))[0]


colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

plt.figure(figsize=(13, 9.5))
plt.subplots_adjust(left=.001, right=.999, bottom=.001, top=.96, wspace=.05,
                    hspace=.01)

plot_num = 1

for i_dataset, dataset in enumerate(datasets):
    X, true_labels = dataset
    # normalize dataset for nicer plotting
    X = StandardScaler().fit_transform(X)

    for name, func_choose_nb_cluster in [
            ('Silhouette', max_silhouette),
            ('Stability', stability),
            ('Gap statistic', gap_statistic),
            ('Calinski-Harabasz index', max_CH_index),
            ('Distortion jump', distortion_jump),
    ]:
        # predict cluster memberships
        t0 = time.time()
        nb_cluster = func_choose_nb_cluster(X, cluster_estimator, k_max=10)
        t1 = time.time()

        # retrieving clustering done
        cluster_estimator.set_params(n_clusters=nb_cluster)
        y_pred = cluster_estimator.fit_predict(X)

        # plot
        plt.subplot(3, 5, plot_num)
        if i_dataset == 0:
            plt.title(name, size=18)
        plt.scatter(X[:, 0], X[:, 1], color=colors[y_pred].tolist(), s=10)

        plt.xlim(-2, 2)
        plt.ylim(-2, 2)
        plt.xticks(())
        plt.yticks(())
        plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
                 transform=plt.gca().transAxes, size=15,
                 horizontalalignment='right')
        plot_num += 1

plt.show()
