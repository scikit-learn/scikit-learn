"""
=================================
Handwritten digits decompositions
=================================

This example compares different unsupervised matrix decomposition (dimension
reduction) methods on the digits dataset.

"""
print __doc__

# Authors: Vlad Niculae, Alexandre Gramfort
# License: BSD

from time import time

import pylab as pl

from scikits.learn.decomposition import RandomizedPCA, NMF, SparsePCA, FastICA
from scikits.learn.cluster import KMeans
from scikits.learn.datasets import load_digits

n_row, n_col = 4, 4
n_components = n_row * n_col

###############################################################################
# Load digits data
digits = load_digits()
threes = digits.data[digits.target == 3]
threes_centered = threes - threes.mean(axis=0)
print "Dataset consists of %d images" % len(threes)


###############################################################################
def plot_digit_gallery(title, images):
    pl.figure(figsize=(1. * n_col, 1.13 * n_row))
    pl.suptitle(title, size=16)
    vmax = max(images.max(), -images.min())
    for i, comp in enumerate(images):
        pl.subplot(n_row, n_col, i + 1)
        pl.imshow(comp.reshape((8, 8)), cmap=pl.cm.BrBG,
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
                                          init='nndsvd', beta=5, tol=1e-2,
                                          sparseness='components'),
                                      False, False, False),
    ('independent components (ICA)', FastICA(n_components=n_components,
                                             whiten=True),
                                     True, True, False),
    ('sparse components (SparsePCA)', SparsePCA(n_components=n_components,
                                                alpha=5, tol=1e-4),
                                      True, False, False),
    ('cluster centers (KMeans)', KMeans(k=n_components),
                                 True, False, True)
    ]

###############################################################################
# Do the estimation and plot it
for name, estimator, center, transpose, cluster in estimators:
    print "Extracting the top %d %s..." % (n_components, name)
    t0 = time()
    data = threes
    if center:
        data = threes_centered
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
