import numpy as np
from scipy.misc import imread

from sklearn.feature_extraction.image import extract_patches_2d
from sklearn.decomposition import PCA
from sklearn import datasets

from joblib import Parallel, delayed
import multiprocessing

from datetime import datetime

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
startTime = datetime.now()

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

print 'It took', time.time()-start, 'seconds.'
