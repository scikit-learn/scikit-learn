# -*- coding: utf-8 -*-
"""
============================
Demo of MSTCluster algorithm
============================

Creates a minimum spanning tree of the data, then cuts edges with a weight
higher than a given threshold. The remaining connected components form the
clusters of the data

This example shows this for a small dataset. The minimum spanning tree is
computed and shown using edges connecting the points. Green edges are under
the threshold and therefore are kept. Red edges are over the thrsehold and are
cut. The colours of the remaining points are the final clusters found by the
algorithm.

"""
print(__doc__)

import numpy as np
import pylab as pl
from sklearn.cluster import EAC, MSTCluster
from sklearn.utils.sparsetools import minimum_spanning_tree
from sklearn.metrics import pairwise_distances
from sklearn import datasets

threshold = 5.0
n_samples = 50
X, y = datasets.make_blobs(n_samples=n_samples, random_state=8)
D = pairwise_distances(X, metric='euclidean')

# Need to compute the spanning tree from the original data, not the model
# This is due to the model automatically removing those edges.
span_tree = minimum_spanning_tree(D)
rows, cols = span_tree.nonzero()


# We could set metric to be 'precomputed' here, and fit(D) on the next line.
# This example shows how to set the metric parameter and fit on X instead.
model = MSTCluster(threshold=threshold, metric='euclidean')
labels = model.fit(X).labels_

fig = pl.figure()
fig.suptitle('Minimum Spanning Tree Clustering')


ax = fig.add_subplot(111)
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

for sample1, sample2 in zip(rows, cols):
    p = zip(X[sample1], X[sample2])
    if D[sample1,sample2] < threshold:
        color = 'g'
    else:
        color = 'r'
    ax.plot(p[0], p[1], color=color, linewidth=2, zorder=1)


ax.scatter(X[:, 0], X[:,1], color=colors[labels].tolist(), s=80, zorder=2)

pl.show()

