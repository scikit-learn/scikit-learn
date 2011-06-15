"""
=========================================
 S-curve example with various LLE methods
=========================================

An illustration of dimensionality reduction
with locally linear embedding and its variants
"""

# Author: Jake Vanderplas -- <vanderplas@astro.washington.edu>

print __doc__

from time import time

import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D

from scikits.learn import manifold, datasets

X, color = datasets.samples_generator.s_curve(1000)
#X, color = datasets.samples_generator.swiss_roll(1000)

n_neighbors = 8
out_dim = 2

methods = ['standard', 'ltsa', 'hessian', 'modified']

fig = pylab.figure(figsize=(8, 12))

try:
    # compatibility matplotlib < 1.0
    ax = fig.add_axes((0.25,0.66,0.4,0.3), projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color)
    ax.view_init(4, -72)
except:
    ax = fig.add_axes((0.25,0.66,0.5,0.3))
    ax.scatter(X[:, 0], X[:, 2], c=color)

ax.set_title('Original Data')

for i, method in enumerate(methods):
    t0 = time()
    Y, err = manifold.locally_linear_embedding(X, n_neighbors, out_dim,
                                               eigen_solver='arpack',
                                               method=method)
    t1 = time()
    print "%s: %.2g sec" % (methods[i], t1 - t0)
    print ' err = %.2e' % err

    ax = fig.add_subplot(323 + i)
    ax.scatter(Y[:, 0], Y[:, 1], c=color)
    ax.set_title("method = %s" % methods[i])

pylab.show()
