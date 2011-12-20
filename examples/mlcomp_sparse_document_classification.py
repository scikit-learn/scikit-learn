"""
========================================================
Classification of text documents: using a MLComp dataset
========================================================

This is an example showing how the scikit-learn can be used to classify
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

The dataset used in this example is the 20 newsgroups dataset and should be
downloaded from the http://mlcomp.org (free registration required):

  http://mlcomp.org/datasets/379

Once downloaded unzip the archive somewhere on your filesystem.
For instance in::

  % mkdir -p ~/data/mlcomp
  % cd  ~/data/mlcomp
  % unzip /path/to/dataset-379-20news-18828_XXXXX.zip

You should get a folder ``~/data/mlcomp/379`` with a file named ``metadata``
and subfolders ``raw``, ``train`` and ``test`` holding the text documents
organized by newsgroups.

Then set the ``MLCOMP_DATASETS_HOME`` environment variable pointing to
the root folder holding the uncompressed archive::

  % export MLCOMP_DATASETS_HOME="~/data/mlcomp"

Then you are ready to run this example using your favorite python shell::

  % ipython examples/mlcomp_sparse_document_classification.py

"""
print __doc__

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

from time import time
import sys
import os
import numpy as np
import scipy.sparse as sp
import pylab as pl

from sklearn.datasets import load_mlcomp
from sklearn.feature_extraction.text import Vectorizer
from sklearn.linear_model.sparse import SGDClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.naive_bayes import MultinomialNB


if 'MLCOMP_DATASETS_HOME' not in os.environ:
    print "MLCOMP_DATASETS_HOME not set; please follow the above instructions"
    sys.exit(0)

# Load the training set
print "Loading 20 newsgroups training set... "
news_train = load_mlcomp('20news-18828', 'train')
print news_train.DESCR
print "%d documents" % len(news_train.filenames)
print "%d categories" % len(news_train.target_names)

print "Extracting features from the dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer()
X_train = vectorizer.fit_transform((open(f).read()
                                    for f in news_train.filenames))
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_train.shape
assert sp.issparse(X_train)
y_train = news_train.target

print "Loading 20 newsgroups test set... "
news_test = load_mlcomp('20news-18828', 'test')
t0 = time()
print "done in %fs" % (time() - t0)

print "Predicting the labels of the test set..."
print "%d documents" % len(news_test.filenames)
print "%d categories" % len(news_test.target_names)

print "Extracting features from the dataset using the same vectorizer"
t0 = time()
X_test = vectorizer.transform((open(f).read() for f in news_test.filenames))
y_test = news_test.target
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_test.shape


###############################################################################
# Benchmark classifiers
def benchmark(clf_class, params, name):
    print "parameters:", params
    t0 = time()
    clf = clf_class(**params).fit(X_train, y_train)
    print "done in %fs" % (time() - t0)

    if hasattr(clf, 'coef_'):
        print "Percentage of non zeros coef: %f" % (
          np.mean(clf.coef_ != 0) * 100)

    print "Predicting the outcomes of the testing set"
    t0 = time()
    pred = clf.predict(X_test)
    print "done in %fs" % (time() - t0)

    print "Classification report on test set for classifier:"
    print clf
    print
    print classification_report(y_test, pred,
       target_names=news_test.target_names)

    cm = confusion_matrix(y_test, pred)
    print "Confusion matrix:"
    print cm

    # Show confusion matrix
    pl.matshow(cm)
    pl.title('Confusion matrix of the %s classifier' % name)
    pl.colorbar()


print "Testbenching a linear classifier..."
parameters = {
    'loss': 'hinge',
    'penalty': 'l2',
    'n_iter': 50,
    'alpha': 0.00001,
    'fit_intercept': True,
}

benchmark(SGDClassifier, parameters, 'SGD')

print "Testbenching a MultinomialNB classifier..."
parameters = {'alpha': 0.01}

benchmark(MultinomialNB, parameters, 'MultinomialNB')

pl.show()
