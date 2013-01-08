"""
Online learning of a dictionary of parts of faces
==================================================

This example uses a large dataset of faces to learn a set of 20 x 20
images patches that constitute faces.

From the programming standpoint, it is interesting because it shows how
to use the online API of the scikit-learn to process a very large
dataset by chunks. The way we proceed is that we load an image at a time
and extract randomly 15 patches from this image. Once we have accumulated
750 of these patches (using 50 images), we run the `partial_fit` method
of the online KMeans object, MiniBatchKMeans.
"""
print __doc__

import time

import pylab as pl
import numpy as np


from sklearn import datasets
from sklearn.cluster import MiniBatchKMeans
from sklearn.feature_extraction.image import extract_patches_2d

faces = datasets.fetch_lfw_people()
data = faces.data

###############################################################################
# Learn the dictionary of images

print 'Learning the dictionary... '
rng = np.random.RandomState(0)
kmeans = MiniBatchKMeans(n_clusters=100)
patch_size = (20, 20)

buffer = []
index = 1
t0 = time.time()

# The online learning part
for index, img in enumerate(faces.images):
    data = extract_patches_2d(img, patch_size,
                                max_patches=15, random_state=rng)
    data = np.reshape(data, (len(data), -1))
    buffer.append(data)
    index += 1
    if index % 50 == 0:
        data = np.concatenate(buffer, axis=0)
        data -= np.mean(data, axis=0)
        data /= np.std(data, axis=0)
        kmeans.partial_fit(data)
        buffer = []
    if index % 500 == 0:
        print 'Partial fit of %4i out of %i' % (index,
                                                len(faces.images))

dt = time.time() - t0
print 'done in %.2fs.' % dt

###############################################################################
# Plot the results
pl.figure(figsize=(4.2, 4))
for i, patch in enumerate(kmeans.cluster_centers_):
    pl.subplot(10, 10, i + 1)
    pl.imshow(patch.reshape(patch_size), cmap=pl.cm.gray,
              interpolation='nearest')
    pl.xticks(())
    pl.yticks(())

pl.suptitle('Patches of faces\nTrain time %.1fs on %d patches' %
            (dt, len(faces.images)), fontsize=16)
pl.subplots_adjust(0.08, 0.02, 0.92, 0.85, 0.08, 0.23)

pl.show()
