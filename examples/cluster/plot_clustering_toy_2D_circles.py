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
    (8, 4, 35, 300),
)
noise_level = 0.02
rng = np.random.RandomState(42)
circles = []

for (c_x, c_y, r, n_points) in circle_parameters:
    t = rng.uniform(12 * np.pi, size=n_points)
    noise_x = rng.normal(loc=0.0, scale=noise_level * r, size=n_points)
    noise_y = rng.normal(loc=0.0, scale=noise_level * r, size=n_points)

    circle_x = r * np.cos(t) + c_x + noise_x
    circle_y = r * np.sin(t) + c_y + noise_y
    circles.append(np.array([circle_x, circle_y]).T)

X = np.concatenate(circles)

# Shuffle the samples to ensure that the algo has no way of cheating
indices = np.arange(X.shape[0])
rng.shuffle(indices)
X = X[indices]

# Plot the raw dataset

pl.figure()
pl.scatter(X[:, 0], X[:, 1])
pl.title("Original dataset")

# Build a heat kernel of euclidean distances as affinity matrix
#distances = euclidean_distances(X, X)
#affinity = np.exp(- (distances ** 2) / (0.1 * distances.mean() ** 2))

# Build a knn graph as affinity matrix
affinity = kneighbors_graph(X, n_neighbors=100).toarray()

labels = spectral_clustering(affinity, k=3, mode='amg')

pl.figure()
for l, c in zip(np.unique(labels), 'rgbcmyk'):
    X_l = X[labels == l, :]
    pl.scatter(X_l[:, 0], X_l[:, 1], color=c)
pl.title("Data clustered by spectral clustering")
pl.show()
