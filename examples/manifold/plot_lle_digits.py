"""
===============================================
Handwritten digits and Locally Linear Embedding
===============================================

An illustration of locally linear embedding on the digits dataset.
"""

# Author: Fabian Pedregosa -- <fabian.pedregosa@inria.fr>
# License: BSD, (C) INRIA 2011

print __doc__

import numpy as np
import pylab as pl
from matplotlib import offsetbox

#----------------------------------------------------------------------
# Locally linear embedding of the digits dataset

from scikits.learn import manifold, datasets
digits = datasets.load_digits(n_class=6)

print "Computing LLE embedding"
X_r, err = manifold.locally_linear_embedding(digits.data, 30, 2, reg=1e-2)
print "Done. Reconstruction error: %g" % err


#----------------------------------------------------------------------
# Scale and visualize the embedding vectors

x_min, x_max = np.min(X_r, 0), np.max(X_r, 0)
X_r = (X_r - x_min) / (x_max - x_min)

ax = pl.subplot(111)
for i in range(digits.data.shape[0]):
    pl.text(X_r[i, 0], X_r[i, 1], str(digits.target[i]),
            color=pl.cm.Set1(digits.target[i] / 10.),
            fontdict={'weight': 'bold', 'size': 9})

if hasattr(offsetbox, 'AnnotationBbox'):
    # only print thumbnails with matplotlib > 1.0
    shown_images = np.array([[1., 1.]])  # just something big
    for i in range(digits.data.shape[0]):
        dist = np.sum((X_r[i] - shown_images) ** 2, 1)
        if np.min(dist) < 4e-3:
            # don't show points that are too close
            continue
        shown_images = np.r_[shown_images, [X_r[i]]]
        imagebox = offsetbox.AnnotationBbox(
            offsetbox.OffsetImage(digits.images[i], cmap=pl.cm.gray_r),
            X_r[i])
        ax.add_artist(imagebox)

pl.xticks([]), pl.yticks([])
pl.show()
