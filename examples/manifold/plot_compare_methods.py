r"""
=========================================
 Comparison of Manifold Learning methods
=========================================

An illustration of dimensionality reduction on the S-curve dataset
with various manifold learning methods.  The methods are as follows:

* LLE : Standard Locally Linear Embedding.
  :func:`scikits.learn.manifold.locally_linear`, ``method = 'standard'``
* LTSA : Local Tangent Space Alignment.
  :func:`scikits.learn.manifold.locally_linear`, ``method = 'ltsa'``
* Hessian LLE : Hessian Eigenmapping.
  :func:`scikits.learn.manifold.locally_linear`, ``method = 'hessian``
* Modified LLE : Modified Locally Linear Embedding with multiple weights.
  :func:`scikits.learn.manifold.locally_linear`, ``method = 'modified'``
* Isomap : Isometric Mapping.
  :func:`scikits.learn.manifold.isomap`

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
n_neighbors = 30
out_dim = 2

X, color = datasets.samples_generator.s_curve(n_points)

methods = ['standard', 'ltsa', 'hessian', 'modified']
labels = ['LLE', 'LTSA', 'Hessian LLE', 'Modified LLE']

fig = pylab.figure(figsize=(8, 12))
pylab.suptitle("Manifold Learning with %i points, %i neighbors"
               % (n_points, n_neighbors), fontsize=14)

try:
    # compatibility matplotlib < 1.0
    ax = fig.add_subplot(321, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=pylab.cm.Spectral)
    ax.view_init(4, -72)
except:
    ax = fig.add_subplot(321, projection='3d')
    ax.scatter(X[:, 0], X[:, 2], c=color, cmap=pylab.cm.Spectral)

ax.set_title('Original Data')

for i, method in enumerate(methods):
    t0 = time()
    Y, err = manifold.locally_linear_embedding(
        X, n_neighbors, out_dim, eigen_solver='arpack', method=method)
    t1 = time()
    print "%s: %.2g sec" % (methods[i], t1 - t0)
    print ' err = %.2e' % err

    ax = fig.add_subplot(322 + i)
    ax.scatter(Y[:, 0], Y[:, 1], c=color, cmap=pylab.cm.Spectral)
    ax.set_title("%s (%.2g sec)" % (labels[i], t1 - t0))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())

t0 = time()
Y = manifold.isomap(X, n_neighbors, out_dim)
t1 = time()
print "Isomap: %.2g sec" % (t1 - t0)
ax = fig.add_subplot(326)
ax.scatter(Y[:, 0], Y[:, 1], c=color, cmap=pylab.cm.Spectral)
ax.set_title("Isomap (%.2g sec)" % (t1-t0))
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())

pylab.show()
