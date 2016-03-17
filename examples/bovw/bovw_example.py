import numpy as np
from scipy.misc import imread

from sklearn.feature_extraction.image import extract_patches_2d
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

from joblib import Parallel, delayed

import time

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

# Script starts here
start = time.time()

# Define the parameters in use afterwards
patch_size = (9, 9)
max_patches = 100
n_jobs = -1
n_components = 9
max_patches_classify = 20000
nb_words = 50
rng = np.random.RandomState(0)

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
print 'Extracted patch to build dictionary'
# Create a plain matrix to apply the PCA decomposition
patch_arr = np.array(patch_arr, copy=False)
patch_arr = patch_arr.reshape((patch_arr.shape[0] * patch_arr.shape[1],
                               patch_arr.shape[2]))
# Build a PCA dictionary
dict_PCA = PCA(n_components=n_components)
dict_PCA.fit(patch_arr)
print 'Built the PCA dictionary'

# Extract the data and project them using the PCA codebook
# Extract and project all the image feature
patch_arr = Parallel(n_jobs=n_jobs)(delayed(image_extraction_projection)
                                    (path_im,
                                     dict_PCA,
                                     rng,
                                     patch_size,
                                     max_patches_classify)
                                    for path_im in png_files)
print 'Extracted and projected patches for image classification'
# Apply a stratified K-fold classification in which we will learn
# a dictionary
skf = StratifiedKFold(n_folds=5)

# Get the training and testing index from the first fold
train_idx, test_idx = skf.split(patch_arr, labels).next()

# Build the codebook
# Define the number of words to create the codebook
vq = MiniBatchKMeans(n_clusters=nb_words, verbose=1, init='random',
                     batch_size=10 * nb_words, compute_labels=False,
                     reassignment_ratio=0.0, random_state=1, n_init=3)
# Stack the training example
stack_training = np.vstack([patch_arr[t] for t in train_idx])
# Find the centroids
vq.fit(stack_training)
print 'Codebook learnt'

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

print 'Classification performed - the confusion matrix obtained is:'
print confusion_matrix(test_label, pred)
print 'It took', time.time()-start, 'seconds.'
