"""
========================================
Sparse PCA for digits feature extraction
========================================

:ref:`SparsePCA` extracting the first 12 components out of the subset of
the digits dataset containing only the digit 3.

"""
print __doc__

# Authors: Vlad Niculae, Alexandre Gramfort
# License: BSD

import numpy as np
import pylab as pl

from scikits.learn.decomposition import SparsePCA
from scikits.learn.datasets import load_digits

###############################################################################
# Load data and fit the model
rows, cols = 4, 3
digits = load_digits()
threes = digits.data[digits.target == 3]
threes -= threes.mean(axis=0)  # XXX: use preprocessors
model = SparsePCA(n_components=rows * cols, alpha=5)
model.fit(threes)
span = np.max(np.abs(model.components_))

###############################################################################
# Plot sparse components
fig = pl.figure(figsize=(1.2 * cols, 1.4 * rows))
for i, comp in enumerate(model.components_):
    pl.subplot(rows, cols, i + 1)
    pl.imshow(np.reshape(comp, (8, 8)), interpolation='nearest',
               vmin=-span, vmax=span, cmap=pl.cm.PuOr)
    pl.xticks(())
    pl.yticks(())
pl.subplots_adjust(0.01, 0.15, 0.99, 0.99, 0.04, 0.)
cax = fig.add_axes([0.1, 0.06, 0.8, 0.04])
pl.colorbar(cax=cax, orientation='horizontal')
pl.show()
