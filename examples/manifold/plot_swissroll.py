"""
===================================
Swiss Roll reduction with LLE
===================================

An illustration of Swiss Roll reduction
with locally linear embedding
"""

# Author: Fabian Pedregosa -- <fabian.pedregosa@inria.fr>
# License: BSD, (C) INRIA 2011

print __doc__

import pylab as pl
# This import is needed to modify the way figure behaves
from mpl_toolkits.mplot3d import Axes3D

#----------------------------------------------------------------------
# Locally linear embedding of the swiss roll

from scikits.learn import manifold, datasets
X, color = datasets.samples_generator.swiss_roll(n_samples=1500)

print "Computing LLE embedding"
X_r, err = manifold.locally_linear_embedding(X, n_neighbors=12, out_dim=2)
print "Done. Reconstruction error: %g" % err

#----------------------------------------------------------------------
# Plot result

fig = pl.figure()
try:
    # compatibility matplotlib < 1.0
    ax = fig.add_subplot(211, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color)
except:
    ax = fig.add_subplot(211)
    ax.scatter(X[:, 0], X[:, 2], c=color)

ax.set_title("Original data")
ax = fig.add_subplot(212)
ax.scatter(X_r[:,0], X_r[:,1], c=color)
pl.axis('tight')
pl.xticks([]), pl.yticks([])
pl.title('Projected data')
pl.show()
