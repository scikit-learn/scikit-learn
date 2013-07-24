"""
Agglomerative clustering with and without structure
===================================================

* Demonstrates the percolation effect
* Shows the effect of connectivity graph to capture manifold structure
* Demonstrates the interest of having structure for speed reasons
"""

import time
import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph

# Generate sample data
n_samples = 1500
np.random.seed(0)
t = 1.5 * np.pi * (1 + 3 * np.random.rand(1, n_samples))
x = t * np.cos(t)
y = t * np.sin(t)

X = np.concatenate((x, y))
X += .75 * np.random.randn(2, n_samples)
X = X.T

# Create a graph capturing local connectivity. Larger number of neighbors
# will give more homogeneous clusters to the cost of computation
# time. With a very large number of neighbors, the manifold structure of
# the data is no longer respected
knn_graph = kneighbors_graph(X, 20)

for n_clusters in (30, 4):
    plt.figure(figsize=(12, 7))
    for connectivity in (None, knn_graph):
        for index, linkage in enumerate(('average', 'complete', 'ward')):
            plt.subplot(2, 3, (connectivity is None) * 3 + index + 1)
            model = AgglomerativeClustering(linkage=linkage,
                                            connectivity=connectivity,
                                            n_clusters=n_clusters)
            t0 = time.time()
            model.fit(X)
            elapsed_time = time.time() - t0
            plt.scatter(X[:, 0], X[:, 1], c=model.labels_,
                        cmap=plt.cm.spectral)
            plt.title('linkage=%s, connectivity=%r \n(time %.2fs)' % (linkage,
                      connectivity is not None, elapsed_time),
                      fontdict=dict(verticalalignment='top'))
            plt.axis('off')

        plt.subplots_adjust(bottom=0, top=.94, hspace=0, wspace=0,
                            left=0, right=1)
        plt.suptitle('n_cluster=%i' % n_clusters, size=15)


plt.show()

