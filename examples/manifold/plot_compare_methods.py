r"""
=========================================
 Comparison of Manifold Learning methods
=========================================

An illustration of dimensionality reduction on the S-curve dataset
with various manifold learning methods.

For a discussion and comparison of these algorithms, see the
:ref:`manifold module page <manifold>`
"""

# Author: Jake Vanderplas -- <vanderplas@astro.washington.edu>

print __doc__

from time import time

import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from scikits.learn import manifold, datasets

n_points = 1000
X, color = datasets.samples_generator.make_s_curve(n_points)
n_neighbors = 8
out_dim = 2

def scatter_3D(X, title):
    fig = pylab.figure()
    try:
        # compatibility matplotlib < 1.0
        ax = fig.add_axes((0,0,1,1), projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=pylab.cm.Spectral)
        ax.view_init(4, -72)
    except:
        ax = fig.add_subplot(321, projection='3d')
        ax.scatter(X[:, 0], X[:, 2], c=color, cmap=pylab.cm.Spectral)

    ax.set_title(title)

def scatter_2D(X, title):
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.scatter(Y[:, 0], Y[:, 1], c=color, cmap=pylab.cm.Spectral)
    ax.set_title(title)
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())


scatter_3D(X, "Input data: %i points in 3 dimensions" % n_points)

methods = ['standard', 'ltsa', 'hessian', 'modified']
labels = ['LLE', 'LTSA', 'Hessian LLE', 'Modified LLE']    

for i, method in enumerate(methods):
    t0 = time()
    Y, err = manifold.locally_linear_embedding(
        X, n_neighbors, out_dim, eigen_solver='arpack', method=method)
    t1 = time()

    scatter_2D(Y, "%s (%.2g sec)" % (labels[i], t1 - t0))

t0 = time()
Y = manifold.Isomap(n_neighbors, out_dim).fit_transform(X)
t1 = time()
scatter_2D(Y, "Isomap (%.2g sec)" % (t1 - t0))

pylab.show()
