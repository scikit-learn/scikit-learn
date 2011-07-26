"""
===========================================
 Isomap example with s-curve dataset
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

N = 500
n_neighbors = 15
out_dim = 2

#X, color = datasets.samples_generator.swiss_roll(N)
X, color = datasets.samples_generator.s_curve(N)

print "Computing Isomap embedding on %i points" % N

t0 = time()
Y = manifold.isomap(X, n_neighbors, out_dim)
t1 = time()

print " - completed in %.2g sec" % (t1 - t0)

fig = pylab.figure(figsize=(6, 10))

try:
    # compatibility matplotlib < 1.0
    ax = fig.add_subplot(211, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=pylab.cm.Spectral)
    ax.view_init(4, -72)
except:
    ax = fig.add_subplot(211)
    ax.scatter(X[:, 0], X[:, 2], c=color, cmap=pylab.cm.Spectral)

ax.set_title('Original Data')

pylab.subplot(212)
pylab.scatter(Y[:, 0], Y[:, 1], c=color, cmap=pylab.cm.Spectral)
pylab.title("Isomap Embedding")

pylab.show()
