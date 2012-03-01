#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Connectivity-constrained Cluster
=========================================================

This example displays the use of `Ward clustering` using a
connectivity graph on the `Lena` image. The connectivity
parameter of `Ward` defines for each sample the neighbouring
samples following a given structure. This structure is supplied
by *grid_to_graph*, which is a graph of the pixel-to-pixel
connections, given `lena`'s shape.

"""
print __doc__


# Code source: Gael Varoqueux
# Modified for Documentation merge by Jaques Grobler
# License: BSD

import numpy as np
import pylab as pl
import scipy as sp
from sklearn.feature_extraction.image import grid_to_graph
from sklearn import cluster, datasets

lena = sp.lena()

# Downsample the image by a factor of 4
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
X = np.reshape(lena, (-1, 1))

# the structure of the data: pixels connected to their neighbors
# Graph of the pixel-to-pixel connections
connectivity = grid_to_graph(*lena.shape) 

ward = cluster.Ward(n_clusters=30, connectivity=connectivity)
ward.fit(X) 
labels = np.reshape(ward.labels_, lena.shape)
pl.imshow(labels)
pl.show()
