"""
===================================================
Faces recognition example using eigenfaces and SVMs
===================================================

The dataset used in this example is a preprocessed excerpt of the "Labeled Faces
in the Wild", aka LFW_:

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
from collections import defaultdict

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
# Count occurrences of each category

categories = [f.rsplit('_', 1)[0] for f in face_filenames]

counts = defaultdict(lambda: 0)

for cat in categories:
    counts[cat] += 1

################################################################################
# Index category names into integers suitable for scikit-learn

# TODO: factorize this out as a utility function in scikit-learn

class Vocabulary(dict):

    def __getitem__(self, k):
        if k not in self:
            self[k] = len(self)
        return super(Vocabulary, self).__getitem__(k)

    def add(self, k):
        self[k]

vocabulary = Vocabulary()

for cat in counts.iterkeys():
    vocabulary.add(cat)

category_names = dict((v, k) for k, v in vocabulary.iteritems())


################################################################################
# Sub sample the dataset to restrict to the most frequent categories

target = np.asarray([vocabulary[cat] for cat in categories])

top_categories = [(count, vocabulary[cat])
                  for cat, count in counts.iteritems()]
top_categories.sort(reverse=True)

labels = [i for c, i in top_categories[:5]]
kept = set(labels)

mask = np.asarray([i for i, t in enumerate(target) if t in kept])

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
n_components = 100

print "Extracting the top %d eigenfaces" % n_components
pca = PCA(n_comp=n_components).fit(X_train)

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
print classification_report(y_test, y_pred, labels=labels,
                            class_names=[category_names[l] for l in labels])

print confusion_matrix(y_test, y_pred, labels=labels)


################################################################################
# Qualitative evaluation of the predictions using matplotlib

n_row = 3
n_col = 4

for i in range(n_row * n_col):
    pl.subplot(n_row, n_col, i + 1)
    pl.imshow(X_test[i].reshape((64, 64)), cmap=pl.cm.gray_r)
    pl.title('%s' % category_names[y_test[i]])

pl.show()

# TODO: find a way to hide the x and y axis
# TODO: plot the top eigenfaces and the singular values absolute values

