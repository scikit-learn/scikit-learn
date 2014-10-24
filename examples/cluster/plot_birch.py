# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Manoj Kumar <manojkumarsivaraj334@gmail.com>

# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import Birch
from sklearn.datasets.samples_generator import make_blobs
from itertools import cycle

colors = cycle(['#4EACC5', '#FF9C34', '#4E9A06'])

##############################################################################
# Generate sample data
np.random.seed(0)
batch_size = 45
centers = [[1, 1], [-1, -1], [1, -1]]
n_clusters = len(centers)
X, labels_true = make_blobs(n_samples=3000, centers=centers, cluster_std=0.7)

##############################################################################
for threshold in [3.0, 1.0]:
    # Compute clustering with Birch
    birch = Birch(threshold=threshold, branching_factor=8)
    birch.fit(X)

    ##########################################################################
    # Plot result
    labels = birch.labels_
    centroids = birch.centroids_
    n_clusters = np.unique(labels).size
    print "n_clusters : %d" % n_clusters

    plt.figure()
    for this_centroid, k, col in zip(centroids, range(n_clusters), colors):
        mask = labels == k
        plt.plot(X[mask, 0], X[mask, 1], 'w',
                 markerfacecolor=col, marker='.')
        plt.plot(this_centroid[0], this_centroid[1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=10)
    plt.title('Birch')
    plt.show()
