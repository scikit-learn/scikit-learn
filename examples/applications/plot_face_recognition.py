"""
===================================================
Faces recognition example using eigenfaces and SVMs
===================================================

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", aka LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/

Expected results for the top 5 most represented people in the dataset::

                     precision    recall  f1-score   support

      George_W_Bush       0.84      0.88      0.86       129
       Colin_Powell       0.80      0.84      0.82        58
         Tony_Blair       0.66      0.62      0.64        34
    Donald_Rumsfeld       0.87      0.79      0.83        33
  Gerhard_Schroeder       0.75      0.64      0.69        28

        avg / total       0.81      0.81      0.81       282

"""
print __doc__

import os
from gzip import GzipFile

import numpy as np
import pylab as pl

from scikits.learn.metrics import classification_report
from scikits.learn.metrics import confusion_matrix
from scikits.learn.pca import PCA
from scikits.learn.svm import SVC

################################################################################
# Download the data, if not already on disk

url = "http://dl.dropbox.com/u/5743203/data/lfw_preprocessed.tar.gz"
archive_name = "lfw_preprocessed.tar.gz"
folder_name = "lfw_preprocessed"

if not os.path.exists(folder_name):
    if not os.path.exists(archive_name):
        import urllib
        print "Downloading data, please Wait (58.8MB)..."
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

faces_filename = os.path.join(folder_name, "faces.npy.gz")
filenames_filename = os.path.join(folder_name, "face_filenames.txt")

faces = np.load(GzipFile(faces_filename))
face_filenames = [l.strip() for l in file(filenames_filename).readlines()]

# normalize each picture by centering brightness
faces -= faces.mean(axis=1)[:, np.newaxis]


################################################################################
# Index category names into integers suitable for scikit-learn

# Here we do a little dance to convert file names in integer indices
# (class indices in machine learning talk) that are suitable to be used
# as a target for training a classifier. Note the use of an array with
# unique entries to store the relation between class index and name,
# often called a 'Look Up Table' (LUT).
# Also, note the use of 'searchsorted' to convert an array in a set of
# integers given a second array to use as a LUT.
categories = np.array([f.rsplit('_', 1)[0] for f in face_filenames])

# A unique integer per category
category_names = np.unique(categories)

# Turn the categories in their corresponding integer label
target = np.searchsorted(category_names, categories)

# Subsample the dataset to restrict to the most frequent categories
selected_target = np.argsort(np.bincount(target))[-5:]

# If you are using a numpy version >= 1.4, this can be done with 'np.in1d'
mask = np.array([item in selected_target for item in target])

X = faces[mask]
y = target[mask]

n_samples, n_features = X.shape

print "Dataset size:"
print "n_samples: %d" % n_samples
print "n_features: %d" % n_features

split = n_samples * 3 / 4

X_train, X_test = X[:split], X[split:]
y_train, y_test = y[:split], y[split:]

################################################################################
# Compute a PCA (eigenfaces) on the training set
n_components = 200

print "Extracting the top %d eigenfaces" % n_components
pca = PCA(n_comp=n_components, do_fast_svd=True).fit(X_train)

eigenfaces = pca.components_.T.reshape((n_components, 64, 64))

# project the input data on the eigenfaces orthonormal basis
X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)


################################################################################
# Train a SVM classification model

print "Fitting the classifier to the training set"
clf = SVC(C=100).fit(X_train_pca, y_train, class_weight="auto")


################################################################################
# Quantitative evaluation of the model quality on the test set

y_pred = clf.predict(X_test_pca)
print classification_report(y_test, y_pred, labels=selected_target,
                            class_names=category_names[selected_target])

print confusion_matrix(y_test, y_pred, labels=selected_target)


################################################################################
# Qualitative evaluation of the predictions using matplotlib

n_row = 3
n_col = 4

pl.figure(figsize=(2*n_col, 2.3*n_row))
pl.subplots_adjust(bottom=0, left=.01, right=.99, top=.95, hspace=.15)
for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(X_test[i].reshape((64, 64)), cmap=pl.cm.gray)
    pl.title('pred: %s\ntrue: %s' % (category_names[y_pred[i]],
                                     category_names[y_test[i]]), size=12)
    pl.xticks(())
    pl.yticks(())

pl.show()

# TODO: plot the top eigenfaces and the singular values absolute values

