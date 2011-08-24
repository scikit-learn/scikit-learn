"""
=========================================
Image denoising using dictionary learning
=========================================
An example comparing the effect of reconstructing noisy fragments
of Lena using online :ref:`DictionaryLearning` and various transform methods.

The dictionary is fitted on the non-distorted left half of the image, and
subsequently used to reconstruct the right half.

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
lena = sp.lena() / 256.0

# downsample for higher speed
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
lena /= 16.0
height, width = lena.shape

# Distort the right half of the image
print "Distorting image..."
distorted = lena.copy()
distorted[:, height/2:] += 0.075 * np.random.randn(width, height/2)

# Extract all clean patches from the left half of the image
print "Extracting clean patches..."
patch_size = (4, 4)
data = extract_patches_2d(distorted[:, :height/2], patch_size)
data = data.reshape(data.shape[0], -1)
intercept = np.mean(data, 0)
data -= intercept

###############################################################################
# Learn the dictionary from clean patches
t0 = time()
dico = DictionaryLearningOnline(n_atoms=100, alpha=1e-2, n_iter=300,
                                verbose=True, transform_algorithm='omp')
V = dico.fit(data).components_
dt = time() - t0
print dt

pl.figure(figsize=(4.5, 5))
for i, comp in enumerate(V):
    pl.subplot(10, 10, i + 1)
    pl.imshow(comp.reshape(patch_size), cmap=pl.cm.gray_r)
    pl.xticks(())
    pl.yticks(())
pl.suptitle("Dictionary learned from Lena patches\n" +
            "Train time %.1fs on %d patches" % (dt, len(data)),
            fontsize=16)
pl.subplots_adjust(0.02, 0.01, 0.98, 0.88, 0.08, 0.01)

def show_with_diff(image, reference, title):
    pl.figure(figsize=(5, 3.2))
    pl.subplot(1, 2, 1)
    pl.title("Image")
    pl.imshow(image, vmin=0, vmax=1, cmap=pl.cm.gray, interpolation='nearest')
    pl.xticks(())
    pl.yticks(())
    pl.subplot(1, 2, 2)
    difference = image - reference

    pl.title("Difference (norm: %.2f)" % np.sqrt(np.sum(difference ** 2)))
    pl.imshow(difference, vmin=-0.5, vmax=0.5, cmap=pl.cm.PuOr,
              interpolation='nearest')
    pl.xticks(())
    pl.yticks(())
    pl.suptitle(title, size=16)
    pl.subplots_adjust(0.02, 0.03, 0.98, 0.79, 0.02, 0.2)

###############################################################################
# Display the distorted image
show_with_diff(distorted, lena, "Distorted image")

###############################################################################
# Extract noisy patches and reconstruct them using the dictionary
print "Extracting noisy patches..."
data = extract_patches_2d(distorted[:, height/2:], patch_size, random_state=0)
data = data.reshape(data.shape[0], -1) - intercept

transform_algorithms = [
    ('1-Orthogonal Matching Pursuit', 'omp',
     {'n_nonzero_coefs': 1, 'precompute_gram': True}),

    ('2-Orthogonal Matching Pursuit', 'omp',
     {'n_nonzero_coefs': 2, 'precompute_gram': True}),

    ('5-Least-angle regression', 'lars',
     {'max_iter': 5})]

reconstructions = {}
for title, transform_algorithm, fit_params in transform_algorithms:
    print title,
    reconstructions[title] = lena.copy()
    t0 = time()
    dico.transform_algorithm = transform_algorithm
    code = dico.transform(data, **fit_params)
    patches = np.dot(code, V) + intercept
    patches = patches.reshape(len(data), *patch_size)
    reconstructions[title][:, height/2:] = reconstruct_from_patches_2d(patches,
                                                           (width, height / 2))
    dt = time() - t0
    print dt
    show_with_diff(reconstructions[title], lena,
                   title + ' (time: %.1fs)' % dt)

pl.show()
