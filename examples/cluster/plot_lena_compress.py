#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Vector Quantization Example
=========================================================

The classic image processing example, Lena, an 8-bit grayscale
bit-depth, 512 x 512 sized image, is used here to illustrate
how `k`-means is used for vector quantization.

"""
print(__doc__)


# Code source: Gael Varoqueux
# Modified for Documentation merge by Jaques Grobler
# License: BSD 3 clause

import numpy as np
import scipy as sp
import pylab as pl

from sklearn import cluster

n_clusters = 5
np.random.seed(0)

try:
    lena = sp.lena()
except AttributeError:
    # Newer versions of scipy have lena in misc
    from scipy import misc
    lena = misc.lena()
X = lena.reshape((-1, 1))  # We need an (n_sample, n_feature) array
k_means = cluster.KMeans(n_clusters=n_clusters, n_init=4)
k_means.fit(X)
values = k_means.cluster_centers_.squeeze()
labels = k_means.labels_

# create an array from labels and values
lena_compressed = np.choose(labels, values)
lena_compressed.shape = lena.shape

vmin = lena.min()
vmax = lena.max()

# original lena
pl.figure(1, figsize=(3, 2.2))
pl.imshow(lena, cmap=pl.cm.gray, vmin=vmin, vmax=256)

# compressed lena
pl.figure(2, figsize=(3, 2.2))
pl.imshow(lena_compressed, cmap=pl.cm.gray, vmin=vmin, vmax=vmax)

# equal bins lena
regular_values = np.linspace(0, 256, n_clusters + 1)
regular_labels = np.searchsorted(regular_values, lena) - 1
regular_values = .5 * (regular_values[1:] + regular_values[:-1])  # mean
regular_lena = np.choose(regular_labels.ravel(), regular_values)
regular_lena.shape = lena.shape
pl.figure(3, figsize=(3, 2.2))
pl.imshow(regular_lena, cmap=pl.cm.gray, vmin=vmin, vmax=vmax)

# histogram
pl.figure(4, figsize=(3, 2.2))
pl.clf()
pl.axes([.01, .01, .98, .98])
pl.hist(X, bins=256, color='.5', edgecolor='.5')
pl.yticks(())
pl.xticks(regular_values)
values = np.sort(values)
for center_1, center_2 in zip(values[:-1], values[1:]):
    pl.axvline(.5 * (center_1 + center_2), color='b')

for center_1, center_2 in zip(regular_values[:-1], regular_values[1:]):
    pl.axvline(.5 * (center_1 + center_2), color='b', linestyle='--')

pl.show()
