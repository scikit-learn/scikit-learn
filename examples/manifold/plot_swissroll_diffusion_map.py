# -*- coding: utf-8 -*-
"""
===================================
Swiss Roll reduction
===================================

An illustration of Swiss Roll reduction
with the manifold module (Diffusion Map)

"""

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

################################################################################
# import some data to play with

# The Swissroll dataset
from scikits.learn.datasets import samples_generator
X, Y = samples_generator.swissroll(n_samples=2000)

colors = np.hstack((Y / Y.max(axis=0), np.zeros((len(Y), 1))))

print "Computing Diffusion map embedding"
from scikits.learn.manifold import DiffusionMap
embedding = DiffusionMap(n_coords=2, n_neighbors=8)
X_r = embedding.fit(Y).embedding_
pl.figure()
pl.scatter(X_r[:,0], X_r[:,1], c=colors)
pl.title("Diffusion map reduction")

pl.show()
