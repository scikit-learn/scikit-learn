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
X, Y = samples_generator.swissroll(n_samples=2000)

colors = np.hstack((Y / Y.max(axis=0), np.zeros((len(Y), 1))))

fig = pl.figure()
ax = Axes3D(fig)#fig.gca(projection='3d')
ax.scatter(X[:,0], X[:,1], X[:,2], c = colors)
ax.set_title("Original data")

pl.show()
