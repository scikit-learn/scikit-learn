# -*- coding: utf-8 -*-
"""
===================================
Swiss Roll reduction
===================================

An illustration of Swiss Roll reduction
with the manifold module

"""

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

################################################################################
# import some data to play with

# The Swissroll dataset
from scikits.learn.datasets import samples_generator
X, Y = samples_generator.swissroll()

colors = np.hstack((X / X.max(axis=0), np.zeros((len(X), 1))))

fig = pl.figure()
ax = Axes3D(fig)#fig.gca(projection='3d')
ax.scatter(Y[:,0], Y[:,1], Y[:,2], c = colors)
ax.set_title("Original data")

from scikits.learn.manifold import LLE
embedding = LLE(n_coords=2, n_neighbors=8)
X_r = embedding.fit(X).embedding_
pl.figure()
pl.scatter(X_r[:,0], X_r[:,1], c=colors)
pl.title("LLE reduction")

from scikits.learn.manifold import LaplacianEigenmap
embedding = LaplacianEigenmap(n_coords=2, n_neighbors=8)
X_r = embedding.fit(X).embedding_
pl.figure()
pl.scatter(X_r[:,0], X_r[:,1], c=colors)
pl.title("Laplacian Eigenmap reduction")

pl.show()
