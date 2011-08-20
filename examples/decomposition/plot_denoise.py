"""
=========================================
Image denoising using dictionary learning
=========================================
An example comparing the effect of reconstructing noisy fragments
of Lena using a dictionary learned from the clear image.
"""
print __doc__

from time import time

import pylab as pl
import scipy as sp
import numpy as np

from scikits.learn.decomposition import DictionaryLearningOnline
from scikits.learn.feature_extraction.image import extract_patches_2d, \
                                                   reconstruct_from_patches_2d

###############################################################################
# Load Lena image and extract patches
lena = sp.lena()
lenb = lena.copy()

print "Extracting clean patches..."
data = extract_patches_2d(lena, (6, 6), max_patches=int(1e4), random_state=0)
data = data.reshape(data.shape[0], -1)
intercept = np.mean(data, 0)
data -= intercept

###############################################################################
# Learn the dictionary from clean patches
dico = DictionaryLearningOnline(n_atoms=100, alpha=0.01, n_iter=100,
                                verbose=True, transform_algorithm='omp')
V = dico.fit(data).components_

###############################################################################
# Generate noisy data and reconstruct using various methods
print ""  # a line break
print "Distorting image fragments..."
fragments = [(slice(200, 300), slice(200, 300)),
             (slice(200, 300), slice(300, 400)),
             (slice(300, 400), slice(200, 300)),
             (slice(300, 400), slice(300, 400))]

for i, fragment in enumerate(fragments):
    lena[fragment] += 50 * np.random.randn(100, 100)
    img = lena[fragment]
    data = extract_patches_2d(img, (6, 6))
    data = data.reshape(len(data), -1)
    data -= intercept
    print "Reconstructing image fragment %d..." % (i + 1),
    t0 = time()
    if i == 0:
        code = dico.transform(data, n_nonzero_coefs=1, precompute_gram=True)
    elif i == 1:
        code = dico.transform(data, n_nonzero_coefs=2, precompute_gram=True)
    elif i == 2:
        dico.transform_algorithm = "lars"
        code = dico.transform(data, max_iter=5)
    elif i == 3:
        dico.transform_algorithm = "threshold"
        code = dico.transform(data, alpha=1.0)
    dt = time() - t0
    print " done in %.2fs" % dt
    data = np.dot(code, V) + intercept
    data = data.reshape(len(data), 6, 6)

    lenb[fragment] = reconstruct_from_patches_2d(data, (100, 100))
    if i == 3:  # tresholding breaks the range
        lenb[fragment] -= lenb[fragment].min()
        lenb[fragment] = lenb[fragment] / float(lenb.max()) * 256.0

###############################################################################
# Display the results
vmin, vmax = 0, 256
pl.figure(figsize=(9, 5))
pl.subplot(1, 2, 1)
pl.title("Noisy image")
pl.imshow(lena, vmin=vmin, vmax=vmax, cmap=pl.cm.gray, interpolation='nearest')
pl.xticks(())
pl.yticks(())
pl.subplot(1, 2, 2)
pl.title("Reconstructed image")
pl.imshow(lenb, vmin=vmin, vmax=vmax, cmap=pl.cm.gray, interpolation='nearest')
pl.xticks(())
pl.yticks(())
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
params = dict(ha="center", va="center", size=12, bbox=bbox_props)
pl.text(150, 100, "OMP-1", **params)
pl.text(400, 100, "OMP-2", **params)
pl.text(150, 450, "LARS", **params)
pl.text(400, 450, "Threshold", **params)
pl.suptitle("Image denoising with dictionary learning", fontsize=16)
pl.subplots_adjust(0.02, 0.1, 0.98, 0.84, 0.02, 0.2)
pl.show()
