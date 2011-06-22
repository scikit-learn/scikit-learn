from scikits.learn.decomposition.sparse_pca import dict_learning
from scikits.learn.datasets import load_digits
from scikits.learn.feature_extraction.image import extract_patches_2d

import pylab as pl
import scipy as sp

lena = sp.lena()

data = extract_patches_2d(lena, (8, 8), max_patches=10000, seed=0)
data = data.reshape(data.shape[0], 64)

V = dict_learning(data, 36, 10, n_iter=1000, return_code=False, verbose=1)

pl.figure(figsize=(6, 6))
for i, comp in enumerate(V):
    pl.subplot(6, 6, i + 1)
    pl.imshow(comp.reshape((8, 8)), cmap=pl.cm.gray_r)
    pl.xticks(())
    pl.yticks(())
pl.suptitle("Dictionary learned from Lena patches", fontsize=16)
pl.subplots_adjust(0.05, 0.02, 0.95, 0.91, 0.01, 0.04)
pl.show()