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

from scikits.learn.decomposition.dict_learning import DictionaryLearningOnline
from scikits.learn.feature_extraction.image import extract_patches_2d

###############################################################################
# Load Lena image and extract patches
lena = sp.lena()
data = extract_patches_2d(lena, (4, 4), max_patches=int(1e5), seed=0)
data = data.reshape(data.shape[0], 16)
data -= np.mean(data, 0)
data /= np.std(data, 0)

###############################################################################
# Learn dictionary
#V = dict_learning(data, n_atoms=36, alpha=1e-2, n_iter=3000, return_code=False,
#                  verbose=True)

dico = DictionaryLearningOnline(n_atoms=100, alpha=1e-2, verbose=True)
dico = dico.fit(data)

###############################################################################
# Plot dictionary atoms
pl.figure(figsize=(4.5, 5))
for i, comp in enumerate(dico.components_):
    pl.subplot(10, 10, i + 1)
    pl.imshow(comp.reshape((4, 4)), cmap=pl.cm.gray_r)
    pl.xticks(())
    pl.yticks(())
pl.suptitle("Dictionary learned from Lena patches", fontsize=16)
pl.subplots_adjust(0.02, 0.05, 0.98, 0.92, 0.08, 0.01)

pl.show()