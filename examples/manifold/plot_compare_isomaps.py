"""
=========================================
 Comparison of Isomap methods
=========================================

An illustration of dimensionality reduction on the Swiss-roll
data set with both Isomap and Landmark Isomap (L-Isomap) methods.
"""

# Author: Lucas David -- <ld492@drexel.edu>
# License: BSD 3 clause, (C) 2015


from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D

from sklearn import datasets
from sklearn.manifold import Isomap

Axes3D

print(__doc__)

EXECUTIONS = [{'n_samples': n_samples, 'noise': 0.0}
              for n_samples in (1000, 5000, 10000, 30000)]
N_EXECUTIONS = len(EXECUTIONS)


def plot_manifold(X, y, ax, title):
    if X.shape[1] == 3:
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=plt.cm.Spectral)
        ax.view_init(4, -72)
        ax.zaxis.set_major_formatter(NullFormatter())
    else:
        plt.scatter(X[:, 0], X[:, 1], c=y, cmap=plt.cm.Spectral)
    plt.title(title, fontdict={'fontsize': 9})
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    for spine in ('top', 'bottom', 'left', 'right'):
        ax.spines[spine].set_visible(False)
    plt.axis('tight')


np.random.rand(7421)

fig = plt.figure(figsize=(15, 9))
plt.suptitle('ISOMAP and L-ISOMAP Comparison', fontsize=14)

for i, execution in enumerate(EXECUTIONS):
    n_samples, noise = execution['n_samples'], execution['noise']
    print('comparison %i: %i samples' % (i, n_samples))

    X, color = datasets.make_swiss_roll(n_samples, noise=noise,
                                        random_state=3241)

    # Original ISOMAP will only look to a subset,
    # because of its performance issues.
    Y = X[np.random.choice(X.shape[0], size=1000)]

    isomap_start = time()
    transformer = Isomap(n_neighbors=10, n_components=2)
    isomap_Z = transformer.fit(Y).transform(X)
    isomap_end = time()

    l_isomap_start = time()
    transformer = Isomap(n_neighbors=10, n_components=2, landmarks='auto')
    l_isomap_Z = transformer.fit_transform(X)
    l_isomap_end = time()

    print('\tISOMAP took %.6g seconds to fit %i samples.'
          % (isomap_end - isomap_start, 1000))
    print('\tL-ISOMAP took %.6g seconds to fit %s samples.\n'
          % (l_isomap_end - l_isomap_start, n_samples))

    ax = fig.add_subplot(N_EXECUTIONS, 3, i * 3 + 1, projection='3d')
    plot_manifold(X, color, ax, 'X dataset (%s samples)' % n_samples)

    ax = fig.add_subplot(N_EXECUTIONS, 3, i * 3 + 2)
    plot_manifold(isomap_Z, color, ax,
                  'ISOMAP (%.6g sec)' % (isomap_end - isomap_start))

    ax = fig.add_subplot(N_EXECUTIONS, 3, i * 3 + 3)
    plot_manifold(l_isomap_Z, color, ax,
                  'L-ISOMAP (%.6g sec)' % (l_isomap_end - l_isomap_start))

plt.show()
