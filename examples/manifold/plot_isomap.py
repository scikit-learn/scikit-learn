"""
===========================================
 Isomap example with swissroll dataset
===========================================

An illustration of dimensionality reduction
with Isomap.
"""

# Author: Jake Vanderplas == <vanderplas@astro.washington.edu>

print __doc__

from time import time

import numpy as np
import pylab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from scikits.learn import manifold, datasets

N = 200
n_neighbors = 4
out_dim = 2

#X, color = datasets.samples_generator.swiss_roll(N)
X, color = datasets.samples_generator.s_curve(N)

t0 = time()
Y1 = manifold.isomap(X, n_neighbors, out_dim, 'dense')
t1 = time()
Y2 = manifold.isomap(X, n_neighbors, out_dim, 'arpack')
t2 = time()

print "dense : %.2g sec" % (t1 - t0)
print "arpack: %.2g sec" % (t2 - t1)

fig = pylab.figure(figsize=(10, 10))

try:
    # compatibility matplotlib < 1.0
    ax = fig.add_subplot(221, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=pylab.cm.Spectral)
    ax.view_init(4, -72)
except:
    ax = fig.add_subplot(221)
    ax.scatter(X[:, 0], X[:, 2], c=color, cmap=pylab.cm.Spectral)

ax.set_title('Original Data')

pylab.subplot(223)
pylab.scatter(Y1[:, 0], Y1[:, 1], c=color, cmap=pylab.cm.Spectral)
pylab.title("dense eigensolver")

pylab.subplot(224)
pylab.scatter(Y2[:, 0], Y2[:, 1], c=color, cmap=pylab.cm.Spectral)
pylab.title("arpack eigensolver")

pylab.show()
