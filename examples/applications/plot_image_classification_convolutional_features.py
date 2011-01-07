"""
==================================================================
Image classification example using convolutional features and SVMs
==================================================================

The dataset used in this example is the CIFAR-10 dataset:

  http://www.cs.toronto.edu/~kriz/cifar.html

This implementation uses an unsupervised feature extraction scheme
to extract a dictionnary of 400 small (6, 6)-shaped filters to be
convolationally applied to the input images as described in:

  An Analysis of Single-Layer Networks in Unsupervised Feature Learning
  Adam Coates, Honglak Lee and Andrew Ng. In NIPS*2010 Workshop on
  Deep Learning and Unsupervised Feature Learning.

Expected results:

  TODO

"""
print __doc__

import os
import math
import cPickle
from gzip import GzipFile
from time import time

import numpy as np
import pylab as pl

from scikits.learn.grid_search import GridSearchCV
from scikits.learn.metrics import classification_report
from scikits.learn.metrics import confusion_matrix
from scikits.learn.feature_extraction.image import ConvolutionalKMeansEncoder
from scikits.learn.svm import SVC
from scikits.learn.svm import LinearSVC
from scikits.learn.preprocessing import Scaler

################################################################################
# Download the data, if not already on disk

url = "http://www.cs.toronto.edu/~kriz/cifar-10-python.tar.gz"
archive_name = url.rsplit('/', 1)[1]
folder_name = 'cifar-10-batches-py'

if not os.path.exists(folder_name):
    if not os.path.exists(archive_name):
        import urllib
        print "Downloading data, please Wait (163MB)..."
        print url
        opener = urllib.urlopen(url)
        open(archive_name, 'wb').write(opener.read())
        print

    import tarfile
    print "Decompressiong the archive: " + archive_name
    tarfile.open(archive_name, "r:gz").extractall()
    print

################################################################################
# Load dataset in memory

X_train = []
y_train = []

for filename in sorted(os.listdir(folder_name)):
    filepath = os.path.join(folder_name, filename)
    if filename.startswith('data_batch_'):
        dataset = cPickle.load(file(filepath, 'rb'))
        X_train.append(dataset['data'])
        y_train.append(dataset['labels'])
    elif filename == 'test_batch':
        dataset = cPickle.load(file(filepath, 'rb'))
        X_test = np.asarray(dataset['data'], dtype=np.float32)
        y_test = dataset['labels']
    elif filename == 'batch.meta':
        dataset = cPickle.load(file(filepath, 'rb'))
        label_neams = dataset['label_names']

X_train = np.asarray(np.concatenate(X_train), dtype=np.float32)
y_train = np.concatenate(y_train)

#n_samples = X_train.shape[0]

# restrict training size for faster runtime as a demo
n_samples = 1000
X_train = X_train[:n_samples]
y_train = y_train[:n_samples]

# reshape pictures to there natural dimension
X_train = X_train.reshape((X_train.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)
X_test = X_test.reshape((X_test.shape[0], 3, 32, 32)).transpose(0, 2, 3, 1)

# convert to graylevel images for now
X_train = X_train.mean(axis=-1)
X_test = X_test.mean(axis=-1)
#pl.imshow(X_train[0], interpolation='nearest'); pl.show()

# scale dataset
print "scaling images to centered, unit variance vectors"
scaler = Scaler().fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)


################################################################################
# Extract filters

whiten = True # perform whitening or not before kmeans
n_components = 10 # number of singular vectors to keep when whitening

n_centers = 400 # kmeans centers: convolutional filters
patch_size = 6  # size of the side of one filter
max_iter = 5 # kmeans EM iteration

extractor = ConvolutionalKMeansEncoder(
    n_centers=n_centers, patch_size=patch_size, whiten=whiten,
    n_components=n_components, max_iter=max_iter)

print "training convolutional whitened kmeans feature extractor..."
t0 = time()
extractor.fit(X_train)
print "done in %0.3fs" % (time() - t0)

if whiten:
    vr = extractor.pca.explained_variance_ratio_
    print "explained variance ratios for %d kept PCA components:" % vr.shape[0]
    print vr
print "kmeans remaining inertia: %0.3fe6" % (extractor.inertia_ / 1e6)

################################################################################
# Qualitative evaluation of the extracted filters

filters = extractor.kernels_

#from scikits.learn.feature_extraction.image import extract_patches2d
#filters = extract_patches2d(X_train, (32, 32), (6, 6))

n_row = int(math.sqrt(n_centers))
n_col = int(math.sqrt(n_centers))

pl.figure()
for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(filters[i].reshape((patch_size, patch_size)),
              cmap=pl.cm.gray, interpolation="nearest")
    pl.xticks(())
    pl.yticks(())

pl.show()


