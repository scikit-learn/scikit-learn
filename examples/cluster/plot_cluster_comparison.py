"""
=========================================================
Comparing different clustering algorithms on toy datasets
=========================================================

This example aims at showing characteristics of different
clustering algorithms on datasets that are "interesting"
but still in 2D. The last dataset is an example of a 'null'
situation for clustering: the data is homogeneous, and
there is no good clustering.

While these examples give some intuition about the algorithms,
this intuition might not apply to very high dimensional data.

The results could be improved by tweaking the parameters for
each clustering strategy, for instance setting the number of
clusters for the methods that needs this parameter
specified. Note that affinity propagation has a tendency to
create many clusters. Thus in this example its two parameters
(damping and per-point preference) were set to to mitigate this
behavior.
"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.cluster import KMeans
from sklearn.cluster import Ward
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import euclidean_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.datasets import make_circles, make_moons, make_blobs
from sklearn.preprocessing import Scaler

np.random.seed(0)

# Generate datasets
n_samples = 300
noisy_circles = make_circles(n_samples=n_samples, factor=.5, noise=.05)
noisy_moons = make_moons(n_samples=n_samples, noise=.05)
blobs = make_blobs(n_samples=n_samples, random_state=8)
no_structure = np.random.rand(n_samples, 2), None

colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

pl.figure(figsize=(14, 10))
pl.subplots_adjust(left=.001, right=.999, bottom=.01, top=.95, wspace=.05,
        hspace=.01)

plot_num = 1
for i_dataset, dataset in enumerate([noisy_circles, noisy_moons, blobs,
                no_structure]):
    X, y = dataset
    # normalize dataset for easier parameter selection
    X = Scaler().fit_transform(X)

    # estimate bandwidth for mean shift
    bandwidth = estimate_bandwidth(X, quantile=0.3)

    # connectivity matrix for structured Ward
    connectivity = kneighbors_graph(X, n_neighbors=10)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)

    # Compute distances
    distances = euclidean_distances(X)

    # create clustering estimators
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    two_means = KMeans(k=2)
    ward_five = Ward(n_clusters=2, connectivity=connectivity)
    spectral = SpectralClustering(k=2, mode='arpack')
    dbscan = DBSCAN(eps=.3)
    affinity_propagation = AffinityPropagation(damping=.9)

    for algorithm in [two_means, dbscan, spectral, ms, ward_five,
                      affinity_propagation]:
        # predict cluster memberships
        if algorithm == spectral:
            algorithm.fit(connectivity)
        elif algorithm == affinity_propagation:
            # Set a low preference to avoid creating too many
            # clusters
            algorithm.fit(-distances, p=-20*distances.max())
        else:
            algorithm.fit(X)
        y_pred = algorithm.labels_.astype(np.int)

        # plot
        pl.subplot(4, 6, plot_num)
        if i_dataset == 0:
            pl.title(str(algorithm).split('(')[0])
        pl.scatter(X[:, 0], X[:, 1], color=colors[y_pred].tolist())

        if hasattr(algorithm, 'cluster_centers_'):
            centers = algorithm.cluster_centers_
            center_colors = colors[:len(centers)]
            pl.scatter(centers[:, 0], centers[:, 1], s=100, c=center_colors)
        pl.xlim(-2, 2)
        pl.ylim(-2, 2)
        pl.xticks(())
        pl.yticks(())
        plot_num += 1

pl.show()
