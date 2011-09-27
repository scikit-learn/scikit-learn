"""
==================================================
Dictionary learning with K-Means on natural scenes
==================================================
"""

print  __doc__

from time import time
import logging
import pylab as pl
import numpy as np

from sklearn.datasets import load_sample_images
from sklearn.feature_extraction.image import PatchExtractor
from sklearn.decomposition import KMeansCoder
from sklearn.preprocessing import scale


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
images = load_sample_images().images
images = np.asanyarray(images, dtype=np.float)

n_atoms = 100

print "Extracting image patches"
t0 = time()
extr = PatchExtractor(patch_size=(6, 6), max_patches=25000, random_state=0)
patches = extr.transform(images)
patches = patches.reshape(len(patches), -1)
patches = scale(patches, axis=1, copy=False)
print "done in %0.3fs" % (time() - t0)

print "Extracting %d atoms from %d patches" % (n_atoms, len(patches))
t0 = time()
kc1 = KMeansCoder(n_atoms, max_iter=10, verbose=True, local_contrast=False,
				  whiten=False).fit(patches)
print "done in %0.3fs" % (time() - t0)

print "Extracting %d whitened atoms from %d patches" % (n_atoms, len(patches))
t0 = time()
kc2 = KMeansCoder(n_atoms, max_iter=10, verbose=True, local_contrast=False,
				  whiten=True).fit(patches)
print "done in %0.3fs" % (time() - t0)
print kc2.components_[0]

n_row = n_col = int(np.sqrt(n_atoms))

titles = ("without whitening PCA", "with whitening PCA")

for img_index, components in enumerate((kc1.components_, kc2.components_)):
    pl.figure(figsize=(2, 3.5))
    pl.suptitle("Dictionary learned with K-Means\non natural scenes\n" +
                titles[img_index])
    for i, atom in enumerate(components):
        pl.subplot(n_row, n_col, i + 1)
        pl.imshow(atom.reshape((6, 6, 3)), interpolation="nearest")
        pl.xticks(())
        pl.yticks(())
    pl.subplots_adjust(0.02, 0.03, 0.98, 0.79, 0.14, 0.01)
pl.show()
