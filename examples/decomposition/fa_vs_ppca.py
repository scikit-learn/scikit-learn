"""
====================================
Factor Analysis vs probabilistic PCA
====================================

This example demonstrates the main difference between
FA and PPCA -- the structure of the noise.
"""

print __doc__

# Authors: Vlad Niculae, Alexandre Gramfort, Christian Osendorfer
# License: BSD

from numpy.random import RandomState
import pylab as pl

from sklearn.datasets import fetch_olivetti_faces
from sklearn import decomposition

image_shape = (64, 64)
rng = RandomState(0)
n_row, n_col = 2, 3
n_components = n_row * n_col

###############################################################################
# Load faces data
dataset = fetch_olivetti_faces(shuffle=True, random_state=rng)
faces = dataset.data

n_samples, n_features = faces.shape

# global centering
faces_centered = faces - faces.mean(axis=0)

# local centering
faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)

###############################################################################
def plot_gallery(title, images):
    pl.figure(figsize=(2. * n_col, 2.26 * n_row))
    pl.suptitle(title, size=16)
    for i, comp in enumerate(images):
        pl.subplot(n_row, n_col, i + 1)
        vmax = max(comp.max(), -comp.min())
        pl.imshow(comp.reshape(image_shape), cmap=pl.cm.gray,
                  interpolation='nearest',
                  vmin=-vmax, vmax=vmax)
        pl.xticks(())
        pl.yticks(())
    pl.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

###############################################################################

print "Computing Factor Analysis ..."
fa = decomposition.FactorAnalysis(n_components=n_components)
fa.fit(faces_centered)
plot_gallery("FA components", fa.components_[:n_components])

# A plot showing learned variance per pixel
var = fa.noise_variance_
vmax = max(var.max(), -var.min())
pl.figure()
pl.imshow(var.reshape(image_shape), cmap=pl.cm.gray,
        interpolation='nearest', vmin=-vmax, vmax=vmax)
pl.xticks(())
pl.yticks(())

print "Computing pPCA..."
ppca = decomposition.ProbabilisticPCA(n_components=n_components)
ppca.fit(faces_centered)
plot_gallery("pPCA components", ppca.components_[:n_components])

print "Scores -- FA:%f - pPCA:%f"%(fa.loglike_[-1],ppca.score(faces_centered).sum())
pl.show()
