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
from sklearn.datasets import make_circles, make_moons, make_blobs

# Generate datasets
n_samples = 50
circles = make_circles(n_samples=n_samples, factor=.5)
moons = make_moons(n_samples=n_samples)
noisy_moons = make_moons(n_samples=n_samples, noise=.1)
blobs = make_blobs(n_samples=n_samples)

i = 1
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 5)

for dataset in [circles, moons, noisy_moons, blobs]:
    X, y = dataset
    # estimate bandwidth for mean shift
    bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=100)
    #bandwidth = .2
    print(bandwidth)

    # create clustering estimators
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    two_means = KMeans(k=2)
    ten_means = KMeans(k=10)
    #affinity = AffinityPropagation()
    ward_five = Ward(n_clusters=5)
    ward_two = Ward(n_clusters=2)

    for algorithm in [two_means, ten_means, ms, ward_five]:
        # predict cluster memberships
        algorithm.fit(X)
        y_pred = algorithm.labels_
        pl.subplot(4, 4, i)
        pl.title(str(algorithm).split('(')[0])
        pl.scatter(X[:, 0], X[:, 1], color=colors[y_pred])
        pl.xticks(())
        pl.yticks(())

        i += 1

pl.show()

