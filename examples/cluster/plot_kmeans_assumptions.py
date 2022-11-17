"""
====================================
Demonstration of k-means assumptions
====================================

This example is meant to illustrate situations where k-means will produce
unintuitive and possibly undesirable clusters.

- Incorrect number of blobs: in a real setting there is no uniquely defined
  **true** number of clusters. An appropriate number of clusters has to be
  decided from data-based criteria and knowledge of aim.
- Anisotropically distributed blobs: k-means consists of minimizing sample's
  euclidean distances to the centroid of the cluster they are assigned
  to. As a consequence, k-means is more appropriated for clusters that are
  normally distributed with a spherical covariance matrix.
- Unequal variance: k-means is equivalent to taking the maximum likelihood
  estimator for a "mixture" of k gaussian distributions with the same variances
  but with possibly different means.
- Unevenly sized blobs: there is no theoretical result about k-means that states
  that it requires similar cluster sizes to perform well, yet minimizing
  euclidean distances does mean that the more sparse and high-dimensional the
  problem is, the higher is the need to run the algorithm with different
  centroid seeds to ensure a minimal inertia.

"""

# Author: Phil Roth <mr.phil.roth@gmail.com>
# License: BSD 3 clause

# %%
# Data generation
# ---------------

import numpy as np
from sklearn.datasets import make_blobs

n_samples = 1500
random_state = 170
transformation = [[0.60834549, -0.63667341], [-0.40887718, 0.85253229]]

X, y = make_blobs(n_samples=n_samples, random_state=random_state)
X_aniso = np.dot(X, transformation)  # Anisotropic blobs
X_varied, y_varied = make_blobs(
    n_samples=n_samples, cluster_std=[1.0, 2.5, 0.5], random_state=random_state
)  # Unequal variance
X_filtered = np.vstack(
    (X[y == 0][:500], X[y == 1][:100], X[y == 2][:10])
)  # Unevenly sized blobs

# %%
# Plot results
# ------------

import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

common_params = {
    "n_init": "auto",
    "random_state": random_state,
}

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

y_pred = KMeans(n_clusters=2, **common_params).fit_predict(X)
axs[0, 0].scatter(X[:, 0], X[:, 1], c=y_pred)
axs[0, 0].set_title("Incorrect Number of Blobs")

y_pred = KMeans(n_clusters=3, **common_params).fit_predict(X_aniso)
axs[0, 1].scatter(X_aniso[:, 0], X_aniso[:, 1], c=y_pred)
axs[0, 1].set_title("Anisotropically Distributed Blobs")

y_pred = KMeans(n_clusters=3, **common_params).fit_predict(X_varied)
axs[1, 0].scatter(X_varied[:, 0], X_varied[:, 1], c=y_pred)
axs[1, 0].set_title("Unequal Variance")

y_pred = KMeans(n_clusters=3, **common_params).fit_predict(X_filtered)
axs[1, 1].scatter(X_filtered[:, 0], X_filtered[:, 1], c=y_pred)
axs[1, 1].set_title("Unevenly Sized Blobs")

plt.show()

# %%
# For an example on how to find a correct number of blobs, see the example
# :ref:`sphx_glr_auto_examples_cluster_plot_kmeans_silhouette_analysis.py`.
#
# For an example on how other clustering methods deal with anisotropic or
# unequal variance blobs, see the example
# :ref:`sphx_glr_auto_examples_cluster_plot_cluster_comparison.py`.
#
# For more details on how to deal with unevenly sized blobs, see
# :ref:`kmeans_sparse_high_dim`.

# %%
# Possible solution
# -----------------

from sklearn.mixture import GaussianMixture

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

y_pred = KMeans(n_clusters=3, **common_params).fit_predict(X)
axs[0, 0].scatter(X[:, 0], X[:, 1], c=y_pred)
axs[0, 0].set_title("Incorrect Number of Blobs")

y_pred = GaussianMixture(n_components=3, covariance_type="full").fit_predict(X_aniso)
axs[0, 1].scatter(X_aniso[:, 0], X_aniso[:, 1], c=y_pred)
axs[0, 1].set_title("Anisotropically Distributed Blobs")

y_pred = GaussianMixture(n_components=3, covariance_type="full").fit_predict(X_varied)
axs[1, 0].scatter(X_varied[:, 0], X_varied[:, 1], c=y_pred)
axs[1, 0].set_title("Unequal Variance")

y_pred = KMeans(n_clusters=3, n_init=10, random_state=random_state).fit_predict(
    X_filtered
)
axs[1, 1].scatter(X_filtered[:, 0], X_filtered[:, 1], c=y_pred)
axs[1, 1].set_title("Unevenly Sized Blobs")

plt.show()
