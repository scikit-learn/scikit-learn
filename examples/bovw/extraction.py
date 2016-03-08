import numpy as np
from scipy.misc import imread

from sklearn.feature_extraction.image import extract_patches_2d
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn.cross_validation import StratifiedKFold
from sklearn.cluster import MiniBatchKMeans
# from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

from joblib import Parallel, delayed
import multiprocessing

import time

# Function to extract a set of texton from one image
def proc_image(path_image, patch_size=(9,9), max_patches=10000):

    # Parameter for the patch extraction
    rng = np.random.RandomState(0)

    # Read the current image
    im = imread(path_image)
    # Extract patches
    patch = extract_patches_2d(im, patch_size=patch_size,
                               max_patches=max_patches, random_state=rng)
    return patch.reshape((max_patches, np.prod(patch_size) * len(im.shape)))

# Function to extract a set of texton from one image and make a PCA projection
def image_extraction_projection(path_image, dict_PCA, patch_size=(9,9)):

    # Parameter for the patch extraction
    max_patches = None
    rng = np.random.RandomState(0)

    # Read the current image
    im = imread(path_image)
    # Extract patches
    patch = extract_patches_2d(im, patch_size=patch_size,
                               max_patches=max_patches, random_state=rng)

    return dict_PCA.transform(patch.reshape((patch.shape[0],
                                             np.prod(patch_size) *
                                             len(im.shape))))

############### Script starts here ###############
start = time.time()

# Parameters for the script
patch_size = (9, 9)
max_patches = 100
n_jobs = -1
n_components = 9

# Load the data
png_files, labels = datasets.fetch_tu_darmstadt()

############### Dictionary learning through PCA ###############

# Extract the patch
patch_arr = Parallel(n_jobs=n_jobs)(delayed(proc_image)(path_im, patch_size, max_patches) 
                                    for path_im in png_files)

print 'Extracted patch to build dictionary'

# Create a plain matrix to apply the PCA decomposition
patch_arr = np.array(patch_arr)
patch_arr = patch_arr.reshape((patch_arr.shape[0] * patch_arr.shape[1], 
                               patch_arr.shape[2]))

# Build a PCA model
dict_PCA = PCA(n_components=n_components)
dict_PCA.fit(patch_arr)

print 'Built the PCA dictionary'

############### Feature extraction and projection ################

# Extract and project all the image feature
patch_arr = Parallel(n_jobs=n_jobs)(delayed(image_extraction_projection)
                                    (path_im, dict_PCA, patch_size)
                                    for path_im in png_files)

print 'Extracted and projected patches for image classification'

############### Feature extraction and projection ################

# Apply a stratified K-fold classification in which we will learn
# a dictionary
skf = StratifiedKFold(labels, n_folds=5)

# For each iteration
for it, (train_idx, test_idx) in enumerate(skf):

    print 'Cross-validation iteration #{}'.format(it+1)

    ##### Build the codebook #####
    # Define the number of words to create the codebook
    nb_words = 1000
    vq = MiniBatchKMeans(n_clusters=nb_words, verbose=1, init='random',
                         batch_size=10 * nb_words, compute_labels=False,
                         reassignment_ratio=0.0, random_state=1, n_init=3)
    # vq = KMeans(n_clusters=nb_words, verbose=10, n_init=4, n_jobs=-1)
    # Stack the training example
    stack_training = np.vstack([patch_arr[t] for t in train_idx])
    # Find the centroids
    vq.fit(stack_training)

    print 'Codebook learnt'

    ##### Build the training and testing data #####
    train_data = []
    for tr_im in train_idx:
        train_data.append(np.histogram(vq.predict(patch_arr[tr_im]),
                                  bins=range(nb_words+1),
                                  density=True))
    train_data = np.array(train_data)
    train_data = np.vstack(train_data[:, 0])
    train_label = labels[train_idx]

    test_data = []
    for te_im in test_idx:
        test_data.append(np.histogram(vq.predict(patch_arr[te_im]),
                                      bins=range(nb_words+1),
                                      density=True))
    test_data = np.array(test_data)
    test_data = np.vstack(test_data[:, 0])
    test_label = labels[test_idx]

    ##### Time for classification #####
    rf = RandomForestClassifier(n_estimators=100)
    pred = rf.fit(train_data, train_label).predict(test_data)

    print 'Classification performed'
    print confusion_matrix(test_label, pred)

    print 'It took', time.time()-start, 'seconds.'
