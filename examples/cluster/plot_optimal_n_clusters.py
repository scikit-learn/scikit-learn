"""Plot the results of the gap criterium."""

# Authors: Thierry Guillemot <thierry.guillemot.work@gmail.com>

import time
import numpy as np
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans, OptimalNClusterSearch
from sklearn.datasets import make_blobs
from sklearn.metrics import calinski_harabaz_score, fowlkes_mallows_score
from sklearn.metrics import silhouette_score
from sklearn.utils import check_random_state


n_samples, n_features, random_state = 1000, 2, 1
parameters = {'n_clusters': np.arange(1, 7)}

rng = check_random_state(random_state)
datasets = [
    ('3 clusters', make_blobs(n_samples=n_samples, n_features=2,
                              random_state=random_state, centers=3)),
    ('5 clusters', make_blobs(n_samples=n_samples, n_features=2,
                              random_state=random_state, centers=5)),
    ('random', (rng.rand(n_samples, n_features),
                np.zeros(n_samples, dtype=int))),
]

estimator = KMeans(n_init=10, random_state=0)
searchers = [
    ('Silhouette', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters,
        fitting_process='unsupervised', metric=silhouette_score)),
    ('Calinski', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters,
        fitting_process='unsupervised', metric=calinski_harabaz_score)),
    ('Stability', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters, random_state=0,
        fitting_process='stability', metric=fowlkes_mallows_score)),
    ('Distortion jump', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters,
        fitting_process='distortion_jump')),
    ('Gap', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters, random_state=0,
        fitting_process='gap')),
    ('Pham', OptimalNClusterSearch(
        estimator=estimator, parameters=parameters, fitting_process='pham')),
]

color = 'bgrcmyk'
plt.figure(figsize=(13, 9.5))
plt.subplots_adjust(left=.001, right=.999, bottom=.001, top=.96, wspace=.05,
                    hspace=.01)
for k, (data_name, data) in enumerate(datasets):
    X, _ = data
    for l, (search_name, search) in enumerate(searchers):
        t0 = time.time()
        y = search.fit(X).predict(X)
        t1 = time.time()

        colors = np.array([color[k] for k in y])
        plt.subplot(len(datasets), len(searchers),
                    len(searchers) * k + l + 1)
        if k == 0:
            plt.title(search_name, size=18)
        plt.scatter(X[:, 0], X[:, 1], color=colors, alpha=.6)
        plt.xticks(())
        plt.yticks(())
        plt.axis('equal')
        plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
                 transform=plt.gca().transAxes, size=15,
                 horizontalalignment='right')
plt.show()
