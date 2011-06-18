"""
===============================================
Handwritten digits and Locally Linear Embedding
===============================================

An illustration of locally linear embedding on the digits dataset.
"""

# Authors: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD, (C) INRIA 2011

print __doc__

import numpy as np
import pylab as pl
from matplotlib import offsetbox
from scikits.learn.utils.fixes import qr_economic
from scikits.learn import manifold, datasets, decomposition, lda

digits = datasets.load_digits(n_class=6)
X = digits.data
y = digits.target
n_samples, n_features = X.shape

#----------------------------------------------------------------------
# Random 2D projection using a random unitary matrix

print "Computing random projection"
rng = np.random.RandomState(42)
Q, _ = qr_economic(rng.normal(size=(n_features, 2)))
X_projected = np.dot(Q.T, X.T).T


#----------------------------------------------------------------------
# Projection on to the first 2 principal components

print "Computing PCA projection"
X_pca = decomposition.RandomizedPCA(n_components=2).fit_transform(X)

#----------------------------------------------------------------------
# Projection on to the first 2 linear discriminant components

print "Computing LDA projection"
X2 = X.copy()
X2.flat[::X.shape[1] + 1] += 0.01 # Make X invertible
X_lda = lda.LDA(n_components=2).fit_transform(X2, y)


#----------------------------------------------------------------------
# Locally linear embedding of the digits dataset

from scikits.learn import manifold, datasets
digits = datasets.load_digits(n_class=6)
method = 'standard'

print "Computing LLE embedding (method = %s)" % method
X_lle, err = manifold.locally_linear_embedding(digits.data, 30, 2, reg=1e-2,
                                             method=method)
print "Done. Reconstruction error: %g" % err

#----------------------------------------------------------------------
# Scale and visualize the embedding vectors

def plot_embedding(X, title=None):
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)

    pl.figure()
    ax = pl.subplot(111)
    for i in range(digits.data.shape[0]):
        pl.text(X[i, 0], X[i, 1], str(digits.target[i]),
                color=pl.cm.Set1(digits.target[i] / 10.),
                fontdict={'weight': 'bold', 'size': 9})

    if hasattr(offsetbox, 'AnnotationBbox'):
        # only print thumbnails with matplotlib > 1.0
        shown_images = np.array([[1., 1.]])  # just something big
        for i in range(digits.data.shape[0]):
            dist = np.sum((X[i] - shown_images) ** 2, 1)
            if np.min(dist) < 4e-3:
                # don't show points that are too close
                continue
            shown_images = np.r_[shown_images, [X[i]]]
            imagebox = offsetbox.AnnotationBbox(
                offsetbox.OffsetImage(digits.images[i], cmap=pl.cm.gray_r),
                X[i])
            ax.add_artist(imagebox)
    pl.xticks([]), pl.yticks([])
    if title is not None:
        pl.title(title)

plot_embedding(X_projected, "Random Projection of the digits")
plot_embedding(X_pca, "Principal Components projection of the digits")
plot_embedding(X_lda, "Linear Discriminant projection of the digits")
plot_embedding(X_lle, "Locally Linear Embedding of the digits")

pl.show()
