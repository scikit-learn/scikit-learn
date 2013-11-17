"""
Agglomerative clustering with different metrics
===============================================

Demonstrates the effect of different metrics on the hierarchical clustering.
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import AgglomerativeClustering
from sklearn.datasets import make_blobs
from sklearn.metrics import pairwise_distances

np.random.seed(0)

# Generate sample data
n_samples = 300
n_features = 500
t = np.pi * np.linspace(0, 1, n_features)
#centers = 10 * np.array([np.cos(3 * t),
#                         np.sin(t),
#                         np.cos(2 * t),
#                         np.sign(np.cos(4 * t)),
#                         .3 * np.sign(np.cos(4 * t)),
#                         ])
from scipy import fftpack
centers = list()
for i in range(1, 6):
    f = np.zeros_like(t)
    f[2*i] = 1
    centers.append(5 * fftpack.idct(f))

X, y = make_blobs(n_samples=n_samples, n_features=n_features,
                  cluster_std=.5, centers=centers)

# Add non-structured noise
#X_noise = 24 * np.random.rand(int(.02 * n_samples), n_features) - 12
#X_noise = 2 * np.random.rand(int(.01 * n_samples), n_features) - 1

#X_noise = 5 * np.random.randn(int(.02 * n_samples), n_features)
X_noise = np.random.rand(int(.02 * n_samples),
                            20 + n_features)
X_noise = .75 * fftpack.idct(X_noise)
X_noise = X_noise[:, 20:]

X = np.concatenate((X, X_noise))
y = np.concatenate((y, len(centers) * np.ones(len(X_noise))))

n_clusters = 5

for index, metric in enumerate(["euclidean", "cityblock",]):# "cosine"]):
    model = AgglomerativeClustering(n_clusters=n_clusters,
                                    linkage="average", affinity=metric)
    model.fit(X)
    plt.figure()
    plt.axes([0, 0, 1, 1])
    for l, c in zip(np.arange(model.n_clusters), 'rgbmyk'):
        plt.plot(X[model.labels_ == l].T, c=c, alpha=.5)
    plt.axis('off')
    plt.suptitle("affinity=%s" % metric, size=18)

for index, metric in enumerate(["euclidean", "cityblock"]):# "cosine"]):
    avg_dist = np.zeros((len(centers) + 1, len(centers) + 1))
    plt.figure(figsize=(5, 5))
    kwargs = dict()
    for i in range(len(centers) + 1):
        for j in range(len(centers) + 1):
            avg_dist[i, j] = pairwise_distances(X[y == i], X[y == j],
                                                metric=metric,
                                                **kwargs).mean()
    avg_dist /= avg_dist.mean()
    for i in range(len(centers) + 1):
        for j in range(len(centers) + 1):
            plt.text(i, j, '%5.3f' % avg_dist[i, j],
                     verticalalignment='center',
                     horizontalalignment='center')

    plt.imshow(avg_dist, interpolation='nearest', cmap=plt.cm.gnuplot2,
               vmin=0)
    plt.suptitle("Interclass %s distances" % metric, size=18)

plt.show()
