"""
=========================================================
Comparing different clustering algorithms on toy datasets
=========================================================

This example aims at showing characteristics of different
clustering algorithms on datasets that are "interesting"
but still in 2D.

While these examples give some intuition about the algorithms,
this intuition might not apply to very high dimensional data.
"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.cluster import KMeans
from sklearn.cluster import Ward
from sklearn.cluster import SpectralClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.datasets import make_circles, make_moons, make_blobs

# Generate datasets
n_samples = 300
circles = make_circles(n_samples=n_samples, factor=.5)
noisy_circles = make_circles(n_samples=n_samples, factor=.5, noise=.05)
moons = make_moons(n_samples=n_samples)
noisy_moons = make_moons(n_samples=n_samples, noise=.05)
blobs = make_blobs(n_samples=n_samples)

i = 1
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 5)

for dataset in [circles, noisy_circles, moons, noisy_moons, blobs]:
    X, y = dataset
    # estimate bandwidth for mean shift
    bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=100)

    connectivity = kneighbors_graph(X, n_neighbors=10)

    connectivity = 0.5 * (connectivity + connectivity.T) # make affinity symmetric
    dists = euclidean_distances(X)
    # create clustering estimators
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    two_means = KMeans(k=2)
    ten_means = KMeans(k=10)
    ward_five = Ward(n_clusters=2, connectivity=connectivity)
    spectral = SpectralClustering(k=2, mode='arpack')

    for algorithm in [two_means, ten_means, spectral, ms, ward_five]:
        # predict cluster memberships
        if algorithm == spectral:
            algorithm.fit(connectivity)
        else:
            algorithm.fit(X)
        y_pred = algorithm.labels_
        pl.subplot(5, 5, i)
        pl.title(str(algorithm).split('(')[0])
        pl.scatter(X[:, 0], X[:, 1], color=colors[y_pred])
        pl.xticks(())
        pl.yticks(())

        i += 1

pl.show()

