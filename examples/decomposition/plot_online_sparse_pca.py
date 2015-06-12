"""
============================
Faces dataset decompositions
============================

This example applies to :ref:`olivetti_faces` online sparse PCA
(a dictionary learning enforcing dictionary atom sparsity),
from the module :py:mod:`sklearn.decomposition` (see the documentation chapter
:ref:`decompositions`), and display some convergence curve.

"""
print(__doc__)

# Authors: Vlad Niculae, Alexandre Gramfort
# License: BSD 3 clause

import logging
from time import time

from numpy.random import RandomState
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt

from sklearn.datasets import fetch_olivetti_faces
from sklearn.decomposition.dict_learning import sparse_encode,\
    MiniBatchDictionaryLearning

import numpy as np

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
n_row, n_col = 5, 6
n_components = n_row * n_col
image_shape = (64, 64)
rng = RandomState(0)

###############################################################################
# Load faces data
dataset = fetch_olivetti_faces(shuffle=True, random_state=rng)
faces = dataset.data

n_samples, n_features = faces.shape

# global centering
faces_centered = faces - faces.mean(axis=0)

# local centering
faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)

print("Dataset consists of %d faces" % n_samples)

###############################################################################
def plot_gallery(title, images, n_col=n_col, n_row=n_row):
    plt.figure(figsize=(1. * n_col, 1.13 * n_row))
    plt.suptitle(title, size=16)
    for i, comp in enumerate(images):
        plt.subplot(n_row, n_col, i + 1)
        vmax = max(comp.max(), -comp.min())
        plt.imshow(comp.reshape(image_shape), cmap=plt.cm.gray,
                   interpolation='nearest',
                   vmin=-vmax, vmax=vmax)
        plt.xticks(())
        plt.yticks(())
    plt.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

###############################################################################
# It is necessary to add regularisation to sparse encoder (either l1 or l2).
# XXX: This should be mentionned in the documentation
dict_learning = MiniBatchDictionaryLearning(n_components=n_components,
                                            alpha=0.1,
                                            n_iter=400, batch_size=10,
                                            fit_algorithm='lars',
                                            transform_algorithm='lasso_lars',
                                            transform_alpha=0.1,
                                            tol=1e-4,
                                            verbose=10,
                                            l1_gamma=0,
                                            random_state=rng,
                                            n_jobs=3,
                                            debug_info=True)

###############################################################################
# Plot a sample of the input data

plot_gallery("First centered Olivetti faces", faces_centered[:n_components])

plt.savefig('faces.pdf')

###############################################################################
# Do the estimation and plot it
name = "Online Dictionary learning"
print("Extracting the top %d %s..." % (n_components, name))
t0 = time()
data = faces
dict_learning.fit(faces_centered)
train_time = (time() - t0)
print("done in %0.3fs" % train_time)
plot_gallery('%s - Train time %.1fs' % (name, train_time),
             dict_learning.components_[:n_components])

np.save('values', dict_learning.values_)
np.save('density', dict_learning.density_)
np.save('residuals', dict_learning.residuals_, )

print("%s - Component density" % name)
print(1 - np.sum(dict_learning.components_ == 0)\
          / float(np.size(dict_learning.components_)))
print("Code density")
code = dict_learning.transform(faces_centered)
print 1 - np.sum(code == 0) / float(np.size(code))
plot_gallery('%s - Reconstruction' % name,
             code[:n_components].dot(dict_learning.components_))

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6, 9))
ax1.plot(dict_learning.values_)
ax1.set_xlabel("Iteration")
ax1.set_ylabel("Dictionary pixel values")

ax2.plot(dict_learning.residuals_)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Objective surrogate")

ax3.plot(dict_learning.density_)
ax3.set_xlabel("Iteration")
ax3.set_ylabel("Dictionary density")

plt.show()
plt.close()
