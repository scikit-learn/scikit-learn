# -*- coding: utf-8 -*-
"""
===================================
Non-Linear dimensionality reduction
===================================

An illustration of non-linear dimensionality reduction
with the manifold module

"""

import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

##############################################################################
# import some data to play with

# The IRIS dataset
from scikits.learn import datasets
iris = datasets.load_iris()

X = iris.data
y = iris.target
target_names = iris.target_names


def display_embedding(embedding, title):
    fig = pl.figure()
    ax = Axes3D(fig)
    for c, i, m, target_name in zip("rgb", [0, 1, 2], ['o', '^', 's'],
        target_names):
        ax.plot3D(embedding[y==i, 0], embedding[y==i, 1], embedding[y==i, 2],
                  m, c=c, label=target_name)
    ax.legend()
    ax.set_title(title)

##############################################################################
# LLE
print "Computing LLE embedding"
from scikits.learn.manifold import LLE
embedding = LLE(n_coords=3, n_neighbors=28)
X_r = embedding.fit(X).embedding_

display_embedding(X_r, 'LLE embedding of IRIS dataset')

##############################################################################
# Laplacian Eigenmap
print "Computing Laplacian Eigenmap embedding"
from scikits.learn.manifold import LaplacianEigenmap
embedding = LaplacianEigenmap(n_coords=3, n_neighbors=30)
X_r = embedding.fit(X).embedding_

display_embedding(X_r, 'Laplacian Eigenmap embedding of IRIS dataset')

##############################################################################
# Diffusion Map
print "Computing Diffusion Map embedding"
from scikits.learn.manifold import DiffusionMap
embedding = DiffusionMap(n_coords=3, kernel_width=10)
X_r = embedding.fit(X).embedding_

display_embedding(X_r, 'Diffusion Map embedding of IRIS dataset')

##############################################################################
# Hessian Map
print "Computing Hessian Map embedding"
from scikits.learn.manifold import HessianMap
embedding = HessianMap(n_coords=3)
X_r = embedding.fit(X).embedding_

display_embedding(X_r, 'Hessian Map embedding of IRIS dataset')

pl.show()
