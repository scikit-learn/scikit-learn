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

import numpy as np
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

######################################################################
# Compute a PCA (eigendigits) on the digit dataset
print "Extracting the top %d eigendigits..." % n_components,
t0 = time()
pca = RandomizedPCA(n_components=n_components, whiten=True)
pca.fit(threes_centered)
print "done in %0.3fs" % (time() - t0)

eigendigits = pca.components_

######################################################################
# Compute a NMF on the digit dataset
print "Extracting %d non-negative features..." % n_components,
t0 = time()
nmf = NMF(n_components=n_components, init='nndsvd', beta=5, tol=1e-2,
          sparseness="components")
nmf.fit(threes)
print "done in %0.3fs" % (time() - t0)

nmfdigits = nmf.components_

######################################################################
# Compute a ICA on the digit dataset
print "Extracting the top %d independent components..." % n_components,
t0 = time()
ica = FastICA(n_components=n_components, whiten=True)
ica.fit(threes_centered.T)
print "done in %0.3fs" % (time() - t0)

icadigits = ica.components_.T

#####################################################################
# Compute Sparse PCA on the digits
print "Extracting %d sparse components..." % n_components,
t0 = time()
spca = SparsePCA(n_components=n_components, alpha=5, tol=1e-4)
spca.fit(threes_centered)
print "done in %0.3fs" % (time() - t0)

spcadigits = spca.components_
spcaspan = np.max(np.abs(spcadigits))

######################################################################
# Compute a K-Means (cluster centers) on the digit dataset
print "Extracting %d cluster centers..." % n_components,
t0 = time()
km = KMeans(k=n_components)
km.fit(threes_centered)
print "done in %0.3fs" % (time() - t0)

kmeansdigits = km.cluster_centers_


######################################################################
# Plot the results
def plot_digit_gallery(title, images):
    pl.figure(figsize=(1. * n_col, 1.13 * n_row))
    pl.suptitle(title, size=16)
    for i, comp in enumerate(images):
        pl.subplot(n_row, n_col, i + 1)
        pl.imshow(comp.reshape((8, 8)), cmap=pl.cm.gray_r,
                  interpolation='nearest')
        pl.xticks(())
        pl.yticks(())
    pl.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)


plot_digit_gallery('Principal components', eigendigits)
plot_digit_gallery('Independent components', icadigits)
plot_digit_gallery('Non-negative components', nmfdigits)
plot_digit_gallery('K-Means cluter centers', kmeansdigits)

###############################################################################
# Plot sparse components (it's a little bit different)

fig = pl.figure(figsize=(1. * n_col, 1.3 * n_row))
for i, comp in enumerate(spcadigits):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(np.reshape(comp, (8, 8)), interpolation='nearest',
               vmin=-spcaspan, vmax=spcaspan, cmap=pl.cm.PuOr)
    pl.xticks(())
    pl.yticks(())

pl.subplots_adjust(0.01, 0.15, 0.99, 0.92, 0.04, 0.)
cax = fig.add_axes([0.1, 0.06, 0.8, 0.04])
pl.colorbar(cax=cax, orientation='horizontal')
pl.suptitle('Sparse components', size=16)

pl.show()
