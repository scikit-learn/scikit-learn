"""
=================================================
Spectral clustering for non convex cluster shapes
=================================================

"""
print __doc__

# Authors:  Olivier Grisel
# License: BSD

import numpy as np
import pylab as pl

from scikits.learn.cluster import k_means
from scikits.learn.cluster import affinity_propagation
from scikits.learn.cluster import mean_shift
from scikits.learn.cluster import spectral_clustering
from scikits.learn.cluster import Ward
from scikits.learn.cluster import power_iteration_clustering
from scikits.learn.metrics.pairwise import euclidean_distances
from scikits.learn.neighbors import kneighbors_graph

# Generate random samples roughly arranged as nested circles

circle_parameters = (
    # (center_x, center_y, radius, n_points)
    (0, 0, 10, 100),
    (8, 0, 25, 200),
    (8, 4, 55, 300),
)
noise_level = 0.05
rng = np.random.RandomState(42)
circles = []

for (center_x, center_y, radius, n_points) in circle_parameters:
    t = rng.uniform(12 * np.pi, size=n_points)

    circle_x = center_x + radius * np.cos(t)
    circle_y = center_y + radius * np.sin(t)
    circle = np.array([circle_x, circle_y]).T
    noise = rng.normal(scale=noise_level * radius, size=(n_points, 2))

    circles.append(circle + noise)

X = np.concatenate(circles)

# Shuffle the samples to ensure that the algo has no way of cheating
indices = np.arange(X.shape[0])
rng.shuffle(indices)
X = X[indices]

# Utility function to plot the results of the various strategies

def plot_labels(labels, title):
    unique_labels = np.unique(labels)
    for l in unique_labels:
        X_l = X[labels == l, :]
        color = pl.cm.hsv(float(l) / unique_labels.shape[0])
        pl.scatter(X_l[:, 0], X_l[:, 1], color=color)
    pl.title(title)
    pl.xticks(())
    pl.yticks(())

# Plot the raw dataset

pl.figure()
pl.subplot(331)
pl.scatter(X[:, 0], X[:, 1])
pl.title("Original dataset")
pl.xticks(())
pl.yticks(())

# K-Means
_, labels, inertia = k_means(X, k=3)
pl.subplot(332)
plot_labels(labels, "K-Means")
print "K-Means inertia: %f" % inertia

# Mean Shift
_, labels = mean_shift(X, bandwidth=28.0)
pl.subplot(333)
plot_labels(labels, "Mean Shift")

# Build a knn graph as affinity matrix
affinity = kneighbors_graph(X, n_neighbors=10)
affinity = 0.5 * (affinity + affinity.T) # make affinity symmetric

# Affinity propagation
# XXX: I cannot get it to work as expected
#_, labels = affinity_propagation(affinity.toarray(), p=0.5)
#pl.subplot(334)
#plot_labels(labels, "Affinity propagation")

# Ward clustering
labels = Ward(n_clusters=3, connectivity=affinity).fit(X).labels_
pl.subplot(335)
plot_labels(labels, "Ward Clustering")

# Spectral Clustering
# XXX: the spectral clustering results is unstable with the amg-based method
# XXX: we should implement the fast_svd method too
labels = spectral_clustering(affinity, k=3, mode='arpack', rng=rng)
pl.subplot(337)
plot_labels(labels, "Spectral Clustering")

# Power iteration
labels = power_iteration_clustering(
    affinity, k=3, n_vectors=10, tol=1e-5, rng=rng, verbose=False,
    plot_vector=False)
print "Power Iteration Clustering inertia: %f" % inertia
pl.subplot(338)
plot_labels(labels, "Power Iteration")

pl.show()
