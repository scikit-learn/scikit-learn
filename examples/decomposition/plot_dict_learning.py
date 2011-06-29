"""
====================================
Dictionary learning on image patches
====================================

Computes and displays a dictionary learned using :ref:`DictionaryLearning` from
a subset of 4x4 patches extracted from the picture of Lena.

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
data = extract_patches_2d(lena, (4, 4), max_patches=2000, seed=0)
data = data.reshape(data.shape[0], 16)
data -= np.mean(data, 0)
data /= np.std(data, 0)

###############################################################################
# Learn dictionary
#V = dict_learning(data, n_atoms=36, alpha=1e-2, n_iter=3000, return_code=False,
#                  verbose=True)

dico = DictionaryLearning(n_atoms=20, alpha=1e-2, verbose=True, tol=1e-5,
                          max_iter=15)
dico = dico.fit(data)

###############################################################################
# Plot dictionary atoms
pl.figure(figsize=(4.4, 5.5))
for i, comp in enumerate(dico.components_):
    pl.subplot(5, 4, i + 1)
    pl.imshow(comp.reshape((4, 4)), cmap=pl.cm.gray_r)
    pl.xticks(())
    pl.yticks(())
pl.suptitle("Dictionary learned from Lena patches", fontsize=16)
pl.subplots_adjust(0.05, 0.02, 0.95, 0.91, 0.01, 0.04)
pl.show()
