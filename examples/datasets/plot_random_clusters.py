"""
==============================================
Plot randomly generated clustered dataset
==============================================

Plot several randomly generated 2D clustered datasets.
This example illustrates the :func:`datasets.make_blobs`
function.

A total number of samples or a list in which each item is
the number of samples per cluster can be passed as
argument to :func:`datasets.make_blobs`
"""

print(__doc__)

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_blobs

plt.figure(figsize=(8, 8))
plt.subplots_adjust(bottom=.05, top=.9, left=.05, right=.95)

plt.subplot(221)
plt.title("Varying standard deviation and number of samples",
          fontsize='small')
samples = [20, 20, 5]
cluster_stds = np.array([0.05, 0.2, 0.4])
cluster_centers = np.array([[-1.0, -1.0], [0.0, 0.0], [1.0, 1.0]])
X, y = make_blobs(random_state=42, n_samples=samples, n_features=2,
                  centers=cluster_centers, cluster_std=cluster_stds)
plt.scatter(X[:, 0], X[:, 1], marker='o', c=y)

plt.subplot(222)
plt.title("Four clusters with varying number of samples each",
          fontsize='small')
samples = [2, 8, 32, 64]
cluster_centers = [[-5.0, -5.0], [-5.0, 5.0], [5.0, 5.0], [5.0, -5.0]]
X, y = make_blobs(random_state=42, n_samples=samples, n_features=2,
                  centers=cluster_centers, cluster_std=1.0)
plt.scatter(X[:, 0], X[:, 1], marker='o', c=y)

plt.subplot(223)
plt.title("Two clusters with different standard deviations",
          fontsize='small')
samples = [20, 10]
cluster_centers = [[-1.0, -1.0], [1.0, 1.0]]
cluster_stds = [2.0, 0.5]
X, y = make_blobs(random_state=42, n_samples=samples, n_features=2,
                  centers=cluster_centers, cluster_std=cluster_stds)
plt.scatter(X[:, 0], X[:, 1], marker='o', c=y)

plt.subplot(224)
plt.title("Equally distributed clusters",
          fontsize='small')
samples = 60
cluster_centers = 8
X, y = make_blobs(random_state=42, n_samples=samples, n_features=2,
                  centers=cluster_centers, cluster_std=0.9)
plt.scatter(X[:, 0], X[:, 1], marker='o', c=y)

plt.show()
