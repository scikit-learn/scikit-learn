"""
=========================================================
Bisect K-Means and Regular K-Means Performance Comparison
=========================================================

This example shows differences between Regular K-Means algorithm and Bisecting K-Means.

While with increasing n_clusters K-Means re-initializes cluster centers,
Bisecting K-Means shows its divisive nature by splitting clusters obtained from its
previous iterations.

"""
import warnings

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_moons
from sklearn.cluster import BisectKMeans, KMeans
from sklearn.exceptions import EfficiencyWarning
from sklearn.preprocessing import StandardScaler

print(__doc__)

warnings.filterwarnings("ignore", category=EfficiencyWarning)

# Generate sample data
n_samples = 1000
random_state = 0

X, _ = make_moons(n_samples=n_samples,
                  noise=0.4,
                  random_state=random_state)

# Normalize dataset for easier parameter selection
X = StandardScaler().fit_transform(X)

# Number of cluster centers for KMeans and BisectKMeans
n_clusters_list = [2, 3, 4, 5]

# Algorithms to compare
clustering_algorithms = {
    'Bisect K-Means':    BisectKMeans,
    'K-Means':           KMeans,
}

# Colors in RGB format
colors = np.array(['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628'])

# Make subplots for each variant
# fig, axs = plt.subplots(len(n_clusters_list), len(clustering_algorithms))
fig, axs = plt.subplots(len(clustering_algorithms), len(n_clusters_list),
                        figsize=(12, 5))

axs = axs.T

for i, (algorithm_name, algorithm) in enumerate(clustering_algorithms.items()):
    for j, n_clusters in enumerate(n_clusters_list):
        algo = algorithm(n_clusters=n_clusters, random_state=random_state)

        algo.fit(X)

        y_pred = algo.labels_.astype(int)
        centers = algo.cluster_centers_

        axs[j, i].scatter(X[:, 0], X[:, 1], s=10, color=colors[y_pred])
        axs[j, i].set_title(f"{algorithm_name} : {n_clusters} clusters")

        axs[j, i].scatter(centers[:, 0], centers[:, 1], c='#000000', s=20)


# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    ax.set_xticks([])
    ax.set_yticks([])

plt.show()
