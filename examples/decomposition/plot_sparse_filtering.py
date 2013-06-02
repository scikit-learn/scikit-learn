
import numpy as np
import pylab as pl

from sklearn.decomposition import SparseFiltering

from sklearn.feature_extraction.image import extract_patches_2d

from sklearn.datasets import fetch_olivetti_faces

patch_width = 16
n_patches = 25
n_filters = 64
maxfun = 50

###############################################################################
# Load faces data
dataset = fetch_olivetti_faces(shuffle=True)
faces = dataset.data

n_samples, n_features = faces.shape

# global centering
faces_centered = faces - faces.mean(axis=0)

# local centering
faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)
faces_centered = faces_centered.reshape(n_samples, 64, 64)

print("Dataset consists of %d faces" % n_samples)

###############################################################################
patches = [extract_patches_2d(faces_centered[i], (patch_width, patch_width),
                             max_patches=n_patches, random_state=i)
          for i in range(n_samples)]
patches = np.array(patches).reshape(-1, patch_width * patch_width)

###############################################################################
estimator = SparseFiltering(n_filters=n_filters, maxfun=maxfun, iprint=10)
features = estimator.fit_transform(patches)

# Some plotting
pl.figure(0)
for i in range(estimator.w_.shape[0]):
    pl.subplot(int(np.sqrt(n_filters)), int(np.sqrt(n_filters)), i + 1)
    pl.pcolor(estimator.w_[i].reshape(patch_width, patch_width),
              cmap=pl.cm.gray)
    pl.xticks(())
    pl.yticks(())

pl.figure(1)
pl.hist(features.T)
pl.show()
