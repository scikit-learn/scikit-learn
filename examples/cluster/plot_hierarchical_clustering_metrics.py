"""
Agglomerative clustering with different metrics
===============================================

Demonstrates the effect of different metrics on the hierarchical clustering.
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cluster import AgglomerativeClustering

# Generate sample data
n_samples = 1500
np.random.seed(0)
t = 1.5 * np.pi * (1 + 3 * np.random.rand(1, n_samples))
x = t * np.cos(t)
y = t * np.sin(t)

X = np.concatenate((x, y))
X += .75 * np.random.randn(2, n_samples)
X = X.T
X += (8, 8)

n_clusters = 10

for index, metric in enumerate(["euclidean", "cityblock", "cosine"]):
    model = AgglomerativeClustering(n_clusters=n_clusters,
                                    linkage="average", affinity=metric)
    model.fit(X)
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], c=model.labels_,
               cmap=plt.cm.RdYlBu, linewidth=0)
    ax.axis('off')
    ax.set_title("%s" % metric)
