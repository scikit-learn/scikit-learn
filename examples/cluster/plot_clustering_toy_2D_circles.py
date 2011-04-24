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

from scikits.learn.cluster import spectral_clustering
from scikits.learn.metrics.pairwise import euclidean_distances
from scikits.learn.neighbors import kneighbors_graph

# Generate random samples roughly arranged as nested circles

circle_parameters = (
    # (center_x, center_y, radius, n_points)
    (0, 0, 10, 100),
    (8, 0, 25, 300),
    (8, 4, 55, 300),
)
noise_level = 0.02
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

# Plot the raw dataset

pl.figure()
pl.scatter(X[:, 0], X[:, 1])
pl.title("Original dataset")

# Build a knn graph as affinity matrix
affinity = kneighbors_graph(X, n_neighbors=10)
affinity = 0.5 * (affinity + affinity.T) # make affinity symmetric

labels = spectral_clustering(affinity, k=3)

pl.figure()
for l, c in zip(np.unique(labels), 'rgbcmyk'):
    X_l = X[labels == l, :]
    pl.scatter(X_l[:, 0], X_l[:, 1], color=c)
pl.title("Data clustered by spectral clustering")
pl.show()
