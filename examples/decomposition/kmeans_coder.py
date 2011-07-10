"""
====================================================
Dictionary learning with K-Means on faces image data
====================================================

This shows 400 dictionary atoms learned from 6x6 image patches extracted from
the face recognition dataset. The dictionary atoms are learned using
(:ref:`KMeansCoder`), with and respectively without a whitening PCA transform.

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", aka LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/

.. |kc_no_w| image:: /images/plot_kmeans_coder_1.png
    :scale: 50%

.. |kc_w| image:: /images/plot_kmeans_coder_2.png
    :scale: 50%

.. centered:: |kc_no_w| |kc_w|

"""
print __doc__

from time import time
import logging
import pylab as pl

import numpy as np

from scikits.learn.cross_val import StratifiedKFold
from scikits.learn.datasets import fetch_lfw_people
from scikits.learn.feature_extraction.image import PatchExtractor
from scikits.learn.decomposition import KMeansCoder

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


###############################################################################
# Download the data, if not already on disk and load it as numpy arrays

lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)

# reshape the data using the traditional (n_samples, n_features) shape
faces = lfw_people.data
n_samples, h, w = faces.shape

X = faces.reshape((n_samples, h * w))
X -= X.mean(axis=1)[:, np.newaxis]
n_features = X.shape[1]
X = X.reshape((n_samples, h, w))

# the label to predict is the id of the person
y = lfw_people.target
target_names = lfw_people.target_names
n_classes = target_names.shape[0]

print "Total dataset size:"
print "n_samples: %d" % n_samples
print "n_features: %d" % n_features
print "n_classes: %d" % n_classes


###############################################################################
# Split into a training set and a test set using a stratified k fold

train, test = iter(StratifiedKFold(y, k=4)).next()
X_train, X_test = X[train], X[test]
y_train, y_test = y[train], y[test]


###############################################################################
# Compute a PCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset): unsupervised feature extraction / dimensionality reduction
n_atoms = 400

print "Extracting image patches from %d faces" % len(X_train)
t0 = time()
extr = PatchExtractor(patch_size=(6, 6), max_patches=100, random_state=0)
patches = extr.transform(X_train)
print "done in %0.3fs" % (time() - t0)

print "Extracting %d atoms from %d patches" % (
    n_atoms, len(patches))
t0 = time()
kc1 = KMeansCoder(n_atoms, max_iter=5, verbose=True, whiten=False).fit(patches)
print "done in %0.3fs" % (time() - t0)

print "Extracting %d whitened atoms from %d patches" % (
    n_atoms, len(patches))
t0 = time()
kc2 = KMeansCoder(n_atoms, max_iter=5, verbose=True, whiten=True).fit(patches)
print "done in %0.3fs" % (time() - t0)

###############################################################################
# Qualitative evaluation of the extracted filters

n_row = int(np.sqrt(n_atoms))
n_col = int(np.sqrt(n_atoms))
titles = ["without whitening PCA", "with whitening PCA"]

for img_index, components in enumerate((kc1.components_, kc2.components_)):
    pl.figure(figsize=(5, 6))
    pl.suptitle("Dictionary learned with K-Means on the \n LFW dataset " +
                titles[img_index])
    for i, atom in enumerate(components):
        pl.subplot(n_row, n_col, i + 1)
        pl.imshow(atom.reshape((6, 6)), cmap=pl.cm.gray,
                                        interpolation="nearest")
        pl.xticks(())
        pl.yticks(())
    pl.subplots_adjust(0.02, 0.03, 0.98, 0.90, 0.14, 0.01)
    pl.show()
