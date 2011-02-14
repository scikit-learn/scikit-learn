"""
=================================================
Face verification using Labeled Faces in the Wild
=================================================

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", aka LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/

"""

from scikits.learn.datasets import load_lfw_pairs
from scikits.learn.datasets.lfw import scale_face
from scikits.learn.pca import RandomizedPCA

from time import time
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s')

pairs_train = load_lfw_pairs(subset='train')
pairs_test = load_lfw_pairs(subset='test')

n_faces_train, _, h, w, c = pairs_train.data.shape
X_train_0 = pairs_train.data[:, 0, :, :, :].reshape((n_faces_train, -1))
X_train_1 = pairs_train.data[:, 1, :, :, :].reshape((n_faces_train, -1))
y_train = pairs_train.target

n_faces_test, _, h, w, c = pairs_test.data.shape
X_test_0 = pairs_test.data[:, 0, :, :, :].reshape((n_faces_test, -1))
X_test_1 = pairs_test.data[:, 1, :, :, :].reshape((n_faces_test, -1))
y_test = pairs_test.target

# Extract the top eigenfaces over the first half the training set
n_components = 150

print "Extracting the top %d eigenfaces from %d %dx%d samples" % (
    n_components, n_faces_train, w, h)
t0 = time()
pca = RandomizedPCA(n_components=n_components, whiten=True).fit(X_train_0)
print "done in %0.3fs" % (time() - t0)

eigenfaces = pca.components_.T.reshape((n_components, h, w, c))

# project the input data on the eigenfaces orthonormal basis
