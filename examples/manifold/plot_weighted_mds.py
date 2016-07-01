"""
==================================
Weighted Multi-dimensional scaling
==================================

An illustration of the smacof algorithm with weighting. The example considers
the localization of nodes in a wireless network. Because of range limit,
internode distances are known below a threshold, and other distances are
missing. Weighting enables to solve localization, even when data are missing.
"""

# Author: Charles Vanwynsberghe <charles.vanwynsberghe@outlook.fr>
# Licence: BSD

print(__doc__)
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from sklearn.metrics import euclidean_distances

# %% create input data
prng = np.random.RandomState(4)  # set random state

n_node = 15  # Number of nodes

# Node Positions are set randomly
X = np.zeros((n_node, 2))
X[:, 0] = 1*prng.rand(n_node)
X[:, 1] = 3*prng.rand(n_node)

X -= X.mean(axis=0)

# Create distance matrix
D = euclidean_distances(X)

# Add some noise to distance matrix
noise = 0.01 * D * prng.randn(*D.shape)
noise = (noise + noise.T) / 2
Dn = D + noise

# %% Solve node localization problem

# smacof: all pairwise distances are known and used for smacof minimisation
smacof = MDS(n_components=2, n_init=1,
             max_iter=100, dissimilarity='precomputed',
             n_jobs=1, random_state=3)
smacof.fit(Dn)

# weighted smacof: only local connectivity is used,
# ie distances below a threshold
weight = np.ones_like(D)
weight[D > 2] = 0

smacof_w = MDS(n_components=2, n_init=1,
               max_iter=100, dissimilarity='precomputed',
               n_jobs=1, random_state=3)
smacof_w.fit(Dn, weight=weight)

# Align estimations onto ground truth for plot
clf = PCA(n_components=2)
X = clf.fit_transform(X)
smacof.embedding_ = clf.fit_transform(smacof.embedding_)
smacof_w.embedding_ = clf.fit_transform(smacof_w.embedding_)

# %% Plot results
plt.figure(figsize=(15, 10))

# Non-weighted case

# ground truth + internode distances
ax = plt.subplot(221)
ax.set_title("Input data \n\n all pairwise nodes (blue) are known")
ax.axis('equal')
ax.set_xticks([]), ax.set_yticks([])
ax.scatter(*X.T, facecolor='none', s=60)
# Plot the distances
segments = [[X[i, :], X[j, :]]
            for i in range(len(X)) for j in range(len(X))]
lc = LineCollection(segments, colors=([146./255, 197./255, 222./255, 1]))
lc.set_linewidths(0.5 * np.ones(len(segments)))
ax.add_collection(lc)

# ground truth + estimated
plt.subplot(222)
plt.title("Solving localization \n\n By Smacof")
plt.axis('equal')
plt.xticks([]), plt.yticks([])
plt.scatter(*smacof.embedding_.T, marker='+', s=60,
            color='black', label='estimated')
plt.scatter(*X.T, facecolor='none', s=60, label='ground truth')
plt.legend(scatterpoints=1, loc='best', ncol=2)

# Weighted case

# ground truth + internode distances
ax2 = plt.subplot(223)
ax2.set_title("Local distances (blue) are known, \nothers (red) are missing.")
ax2.axis('equal')
ax2.set_xticks([]), ax2.set_yticks([])
ax2.scatter(*X.T, facecolor='none', s=60)
# Plot the known/unknown distances
segments_1 = []
segments_0 = []
for i in range(len(X)):
    for j in range(len(X)):
        if (weight[i, j] == 1):
            segments_1.append([X[i, :], X[j, :]])
        elif (weight[i, j] == 0):
            segments_0.append([X[i, :], X[j, :]])

lc_0 = LineCollection(segments_0, colors=([202./255, 0./255, 32./255, 1]))
lc_1 = LineCollection(segments_1, colors=([146./255, 197./255, 222./255, 1]))
lc_1.set_linewidths(0.5 * np.ones(len(segments_1)))
lc_0.set_linewidths(0.5 * np.ones(len(segments_0)))
ax2.add_collection(lc_0)
ax2.add_collection(lc_1)

# ground truth + estimated
plt.subplot(224)
plt.title("Weighting Smacof enables to solve localization\n with missing data")
plt.axis('equal')
plt.xticks([]), plt.yticks([])
plt.scatter(*smacof_w.embedding_.T, marker='+', s=60,
            color='black', label='estimated')
plt.scatter(*X.T, facecolor='none', s=60, label='ground truth')
plt.legend(scatterpoints=1, loc='best', ncol=2)

plt.show()
