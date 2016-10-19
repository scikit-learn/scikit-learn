"""
===================
Bag of Visual Words
===================

An illustration of a Bag of Visual Words (BoVW) approach using textons for
image recognition [1]_. The proposed solution is a naive and simplified
solution using a limited number of patches and words and do not rely on any
other well known computer vision features. It aims at illustrating the BoVW
under a limited processing time rather than achieving high classification
performance.

References
----------
.. [1] Varma, Manik, and Andrew Zisserman. "A statistical approach to material
       classification using image patch exemplars." IEEE transactions on
       pattern analysis and machine intelligence 31.11 (2009): 2032-2047.

"""
print(__doc__)

# Author: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD

import numpy as np
from scipy.misc import imread

from sklearn.feature_extraction.image import extract_patches_2d
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.externals.joblib import Parallel, delayed

from tudarmstadt import fetch_tu_darmstadt


def proc_image(path_image, rng, patch_size=(9, 9), max_patches=10000):
    """ Function to extract couple of patches from an image """

    # Read the current image
    im = imread(path_image)
    # Extract patches
    patch = extract_patches_2d(im, patch_size=patch_size,
                               max_patches=max_patches, random_state=rng)
    return patch.reshape((max_patches, np.prod(patch_size) * len(im.shape)))


def image_extraction_projection(path_image, dict_PCA, rng,
                                patch_size=(9, 9), max_patches=100000):
    """ Function to extract a couple of patches from an image and apply PCA """

    # Read the current image
    im = imread(path_image)
    # Extract patches
    patch = extract_patches_2d(im, patch_size=patch_size,
                               max_patches=max_patches, random_state=rng)

    return dict_PCA.transform(patch.reshape((patch.shape[0],
                                             np.prod(patch_size) *
                                             len(im.shape))))

# Define the parameters in use afterwards
patch_size = (9, 9)
max_patches = 100
n_jobs = -1
n_components = 9
max_patches_classify = 2000
nb_words = 50
rng = 42

# Load the data
png_files, labels = fetch_tu_darmstadt()

# We will extract some patches from the different images to construct
# a PCA codebook that will be used afterwards
# Extract the patches
patch_arr = Parallel(n_jobs=n_jobs)(delayed(proc_image)(path_im,
                                                        rng,
                                                        patch_size,
                                                        max_patches)
                                    for path_im in png_files)

print('Extracted patch to build dictionary')

# Create a plain matrix to apply the PCA decomposition
patch_arr = np.array(patch_arr, copy=False)
patch_arr = patch_arr.reshape((patch_arr.shape[0] * patch_arr.shape[1],
                               patch_arr.shape[2]))
# Build a PCA dictionary
dict_PCA = PCA(n_components=n_components, random_state=rng)
dict_PCA.fit(patch_arr)

print('Built the PCA dictionary')

# Extract the data and project them using the PCA codebook
# Extract and project all the image feature
patch_arr = Parallel(n_jobs=n_jobs)(delayed(image_extraction_projection)
                                    (path_im,
                                     dict_PCA,
                                     rng,
                                     patch_size,
                                     max_patches_classify)
                                    for path_im in png_files)

print('Extracted and projected patches for image classification')

# Apply a stratified K-fold classification in which we will learn
# a dictionary
skf = StratifiedKFold(n_splits=5, random_state=rng)

# Get the training and testing index from the first fold
train_idx, test_idx = skf.split(patch_arr, labels).next()

# Build the codebook
# Define the number of words to create the codebook
vq = MiniBatchKMeans(n_clusters=nb_words, verbose=False, init='random',
                     batch_size=10 * nb_words, compute_labels=False,
                     reassignment_ratio=0.0, random_state=rng, n_init=3)
# Stack the training example
stack_training = np.vstack([patch_arr[t] for t in train_idx])
# Find the centroids
vq.fit(stack_training)

print('Codebook learnt')

# Build the training and testing data
train_data = []
for tr_im in train_idx:
    train_data.append(np.histogram(vq.predict(patch_arr[tr_im]),
                                   bins=range(nb_words),
                                   density=True))
train_data = np.array(train_data, copy=False)
train_data = np.vstack(train_data[:, 0])
train_label = labels[train_idx]

test_data = []
for te_im in test_idx:
    test_data.append(np.histogram(vq.predict(patch_arr[te_im]),
                                  bins=range(nb_words),
                                  density=True))
test_data = np.array(test_data, copy=False)
test_data = np.vstack(test_data[:, 0])
test_label = labels[test_idx]

# Classification using Random Forest
rf = RandomForestClassifier(n_estimators=10, random_state=rng)
pred = rf.fit(train_data, train_label).predict(test_data)

print('Classification performed - the confusion matrix obtained is:')
print(confusion_matrix(test_label, pred))
