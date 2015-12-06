"""
=========================================
 Comparison of Isomap methods
=========================================

An illustration of dimensionality reduction on the Swiss-roll
data set with both Isomap and Landmark Isomap (L-Isomap) methods.
"""

# Author: Lucas David -- <lucasdavid@comp.ufscar.br>

print(__doc__)

from time import time
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
from sklearn import manifold, datasets

# Next line to silence pyflakes. This import is needed.
Axes3D

isomap_params = {'n_neighbors': 10, 'n_components': 2}
l_isomap_params = {'n_neighbors': 10, 'n_components': 2, 'n_landmarks': 15}

comparisons = [
    (
        n_samples,
        (manifold.Isomap, isomap_params),
        (manifold.Isomap, l_isomap_params),
    ) for n_samples in range(1000, 20000, 5000)]

n_comparisons = len(comparisons)

fig = plt.figure(figsize=(15, 8))
plt.suptitle("Isomap and L-Isomap Comparison", fontsize=14)

for i, comparison in enumerate(comparisons):
    n_samples = comparison[0]
    methods = comparison[1:]

    X, color = datasets.make_swiss_roll(n_samples, random_state=0)

    print('\nComparison %i: %i samples' % (i, n_samples))

    try:
        # compatibility matplotlib < 1.0
        ax = fig.add_subplot(n_comparisons, 3, i * 3 + 1, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=plt.cm.Spectral)
        ax.view_init(4, -72)
    except:
        ax = fig.add_subplot(251, projection='3d')
        plt.scatter(X[:, 0], X[:, 2], c=color, cmap=plt.cm.Spectral)

    plt.title("%i samples." % n_samples)

    for j, (method, params) in enumerate(methods):
        try:
            t0, t1 = time(), None
            Y = method(**params).fit_transform(X)
            t1 = time()

            print("%s(%s): %.6g sec." % (method.__name__,
                                         str(params), (t1 - t0)))

            ax = fig.add_subplot(n_comparisons, 3, i * 3 + 2 + j)
            plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
            plt.title("%s (%.6g sec)" % (method.__name__, (t1 - t0)))

            ax.xaxis.set_major_formatter(NullFormatter())
            ax.yaxis.set_major_formatter(NullFormatter())
            plt.axis('tight')

        except MemoryError as e:
            print("%s: failed." % method.__name__)

plt.show()
