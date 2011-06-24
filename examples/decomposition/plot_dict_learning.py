"""
====================================
Dictionary learning on image patches
====================================

Computes and displays a dictionary learned using :ref:`DictionaryLearning` from
a subset of 8x8 patches extracted from the picture of Lena.

"""
print __doc__

import numpy as np
import pylab as pl
import scipy as sp

from scikits.learn.decomposition import DictionaryLearning
from scikits.learn.decomposition.sparse_pca import dict_learning
from scikits.learn.feature_extraction.image import extract_patches_2d

###############################################################################
# Load Lena image and extract patches
lena = sp.lena()
data = extract_patches_2d(lena, (8, 8), max_patches=500, seed=0)
data = data.reshape(data.shape[0], 64)
data -= np.mean(data, 0)
data /= np.std(data, 0)

###############################################################################
# Learn dictionary
#V = dict_learning(data, n_atoms=36, alpha=1e-2, n_iter=3000, return_code=False,
#                  verbose=True)

dico = DictionaryLearning(n_atoms=12, alpha=1e-2, verbose=True, tol=1e-5,
                          max_iter=15).fit(data)

###############################################################################
# Plot dictionary atoms
pl.figure(figsize=(4.5, 6))
for i, comp in enumerate(dico.components_):
    pl.subplot(4, 3, i + 1)
    pl.imshow(comp.reshape((8, 8)), cmap=pl.cm.gray_r)
    pl.xticks(())
    pl.yticks(())
pl.suptitle("Dictionary learned from Lena patches", fontsize=16)
pl.subplots_adjust(0.05, 0.02, 0.95, 0.91, 0.01, 0.04)
pl.show()
