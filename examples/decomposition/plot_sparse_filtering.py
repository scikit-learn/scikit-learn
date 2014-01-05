"""
===========================================
Sparse filtering on Olivetti faces
===========================================

Unsupervised learning of features for images from the Olivetti faces dataset
using the sparse filtering algorithm. Linear features for sub-patches of the
Olivetti faces are learned using the sparse filtering algorithm. This algorithm
does not try to model the data's distribution but rather to learn features which
are sparsely activated (in the sense that for each image, only a small subset of
features is activated, that each feature is only activated on a small subset of
the examples, and that features are roughly activated equally often). This
sparsity is encoded as an objective function and L-BFGS is used to minimize this
function.

Plotted are the weight matrices of the features (corresponding roughly to gabor
filters) and feature activation histograms.
"""
print(__doc__)

import numpy as np
import pylab as pl

from sklearn.decomposition import SparseFiltering

from sklearn.feature_extraction.image import extract_patches_2d

from sklearn.datasets import fetch_olivetti_faces

patch_width = 16  # Learn features for patches of size patch_width*patch_width
n_patches = 25  # Determines number of random patches extracted from each image
n_features = 64  # How many features are learned
maxfun = 200  # The maximal number of evaluations of the objective function
iprint = 10  # after how many function evaluations is information printed
             # by L-BFGS. -1 for no information

###############################################################################
# Load faces data, normalize faces, and convert 2d structures
dataset = fetch_olivetti_faces(shuffle=True)
faces = dataset.data

n_samples, _ = faces.shape

faces_centered = faces - faces.mean(axis=0)  # global centering

faces_centered -= \
    faces_centered.mean(axis=1).reshape(n_samples, -1)  # local centering

faces_centered = \
    faces_centered.reshape(n_samples, 64, 64)  # Reshaping to 64*64 pixel images

print("Dataset consists of %d faces" % n_samples)

###############################################################################
# Extract n_patches patches randomly from each image
patches = [extract_patches_2d(faces_centered[i], (patch_width, patch_width),
                              max_patches=n_patches, random_state=i)
              for i in range(n_samples)]
patches = np.array(patches).reshape(-1, patch_width * patch_width)

###############################################################################
#
estimator = SparseFiltering(n_features=n_features, maxfun=maxfun, iprint=iprint)
features = estimator.fit_transform(patches)

# #############################################################################
# Plot weights of features
pl.figure(0, figsize=(12, 10))
pl.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.95,
                  wspace=0.1, hspace=0.4)
for i in range(estimator.w_.shape[0]):
    pl.subplot(int(np.sqrt(n_features)), int(np.sqrt(n_features)), i + 1)
    pl.pcolor(estimator.w_[i].reshape(patch_width, patch_width),
              cmap=pl.cm.gray)
    pl.xticks(())
    pl.yticks(())
    pl.title("Feature %4d" % i)

# Plot feature histogram
pl.figure(1)
pl.hist(features.T)
pl.title("Feature activation histogram")

# Plot Lifetime Sparsity histogram
# Lifetime Sparsity: Each feature should only be active for a few examples
pl.figure(2)
activated_features = (features > 0.1).sum(1) / float(features.shape[1])
pl.hist(activated_features)
pl.xlabel("Feature activation ratio over all examples")
pl.title("Lifetime Sparsity Histogram")

# Plot Population Sparsity histogram
# Population Sparsity: Each example should be represented by only a few active
#                      features
pl.figure(3)
activated_features = (features > 0.1).sum(0) / float(features.shape[0])
pl.hist(activated_features)
pl.xlabel("Ratio of active features in example")
pl.title("Population Sparsity Histogram")

pl.show()
