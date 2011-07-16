"""
============================
Faces dataset decompositions
============================

This example compares different unsupervised matrix decomposition (dimension
reduction) methods on the Olivetti faces dataset.

"""
print __doc__

# Authors: Vlad Niculae, Alexandre Gramfort
# License: BSD

import logging
from time import time

import pylab as pl

from scikits.learn.decomposition import RandomizedPCA, NMF, SparsePCA, FastICA
from scikits.learn.cluster import KMeans
from scikits.learn.datasets import fetch_olivetti_faces


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
n_row, n_col = 4, 4
n_components = n_row * n_col
image_shape = (64, 64)

###############################################################################
# Load digits data
dataset = fetch_olivetti_faces()
faces = dataset.data
faces_centered = faces - faces.mean(axis=0)
print "Dataset consists of %d images" % len(faces)


###############################################################################
def plot_digit_gallery(title, images):
    pl.figure(figsize=(1. * n_col, 1.13 * n_row))
    pl.suptitle(title, size=16)
    vmax = max(images.max(), -images.min())
    for i, comp in enumerate(images):
        pl.subplot(n_row, n_col, i + 1)
        pl.imshow(comp.reshape(image_shape), cmap=pl.cm.BrBG,
                  interpolation='nearest',
                  vmin=-vmax, vmax=vmax)
        pl.xticks(())
        pl.yticks(())
    pl.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

###############################################################################
# List of the different estimators, whether to center and transpose the
# problem, and whether the transformer uses the clustering API.
estimators = [
    ('eigendigits (PCA)', RandomizedPCA(n_components=n_components,
                                        whiten=True),
                          True, False, False),
    ('non-negative components (NMF)', NMF(n_components=n_components,
                                          init='nndsvda', beta=1, tol=1e-3,
                                          sparseness='components'),
                                      False, False, False),
    ('independent components (ICA)', FastICA(n_components=n_components,
                                             whiten=True),
                                     True, True, False),
    ('sparse components (SparsePCA)', SparsePCA(n_components=n_components,
                                                alpha=1, tol=1e-4,
                                                verbose=True, max_iter=3),
                                      True, False, False),
    ('cluster centers (KMeans)', KMeans(k=n_components, tol=1e-2, max_iter=3),
                                 True, False, True)
    ]

###############################################################################
# Do the estimation and plot it
for name, estimator, center, transpose, cluster in estimators:
    print "Extracting the top %d %s..." % (n_components, name)
    t0 = time()
    data = faces
    if center:
        data = faces_centered
    if transpose:
        data = data.T
    estimator.fit(data)
    print "done in %0.3fs" % (time() - t0)
    if cluster:
        components_ = estimator.cluster_centers_
    else:
        components_ = estimator.components_
    if transpose:
        components_ = components_.T
    plot_digit_gallery(name, components_)

pl.show()
