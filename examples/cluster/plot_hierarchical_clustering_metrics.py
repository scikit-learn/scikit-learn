"""
Agglomerative clustering with different metrics
===============================================

Demonstrates the effect of different metrics on the hierarchical clustering.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from sklearn.cluster import AgglomerativeClustering
#from sklearn.datasets import make_blobs
from sklearn.metrics import pairwise_distances

np.random.seed(0)

# Generate sample data
n_samples = 300
n_features = 1000
t = np.pi * np.linspace(0, 1, n_features)
#centers = 10 * np.array([np.cos(3 * t),
#                         np.sin(t),
#                         np.cos(2 * t),
#                         np.sign(np.cos(4 * t)),
#                         .3 * np.sign(np.cos(4 * t)),
#                         ])
#from scipy import fftpack

def sqr(x):
    return np.sign(np.cos(x))

centers = list()
X = list()
y = list()
for i, (phi, a) in enumerate([(.5, .15), (.35, .5), (.5, .6), ]):#(.2, .3)]):
    centers.append(a * sqr(2 * np.pi * t))
    for _ in range(100):
        phase_noise = .0 * np.random.normal()
        amplitude_noise = .0 * np.random.normal()
        additional_noise = .0 * np.random.rand()#n_features)
        #additional_noise = ndimage.gaussian_filter1d(additional_noise,
        #                                             sigma=1)
        X.append(12 * ((a + amplitude_noise) * sqr(6 * (t + phi + phase_noise))
                 + additional_noise))
        y.append(i)

X = np.array(X)
y = np.array(y)


n_centers = len(centers)

X_noise = np.random.rand(10, n_features)
X_noise = ndimage.gaussian_filter1d(X_noise, sigma=1, axis=-1)
X_noise /= X_noise.std(axis=-1)[:, np.newaxis]
X_noise *= 0#-5
X_noise -= X_noise.mean(axis=-1)[:, np.newaxis]

X = np.concatenate((X, X_noise))
y = np.concatenate((y, n_centers * np.ones(len(X_noise))))

n_clusters = n_centers

for index, metric in enumerate(["euclidean", ]): #"cityblock", "cosine"]):
    model = AgglomerativeClustering(n_clusters=n_clusters,
                                    linkage="average", affinity=metric)
    model.fit(X)
    plt.figure()
    plt.axes([0, 0, 1, 1])
    for l, c in zip(np.arange(model.n_clusters), 'rgbmyk'):
        plt.plot(X[model.labels_ == l].T, c=c, alpha=.5)
    plt.axis('off')
    plt.suptitle("affinity=%s" % metric, size=18)

for index, metric in enumerate(["euclidean", ]): #"cityblock", "cosine"]):
    avg_dist = np.zeros((n_centers + 1, n_centers + 1))
    plt.figure(figsize=(5, 5))
    kwargs = dict()
    for i in range(n_centers + 1):
        for j in range(n_centers + 1):
            avg_dist[i, j] = pairwise_distances(X[y == i], X[y == j],
                                                metric=metric,
                                                **kwargs).mean()
    avg_dist /= avg_dist.mean()
    for i in range(n_centers + 1):
        for j in range(n_centers + 1):
            plt.text(i, j, '%5.3f' % avg_dist[i, j],
                     verticalalignment='center',
                     horizontalalignment='center')

    plt.imshow(avg_dist, interpolation='nearest', cmap=plt.cm.gnuplot2,
               )#vmin=.84)
    plt.xticks(range(n_centers + 1))
    plt.yticks(range(n_centers + 1))
    plt.colorbar()
    plt.suptitle("Interclass %s distances" % metric, size=18)

plt.show()
