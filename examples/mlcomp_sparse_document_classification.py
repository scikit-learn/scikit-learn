"""
======================================================
Classification of text documents using sparse features
======================================================

This is an example showing how the scikit-learn can be used to classify
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

The dataset used in this example is the 20 newsgroups dataset and should be
downloaded from the http://mlcomp.org (free registration required):

  http://mlcomp.org/datasets/379

Once downloaded unzip the arhive somewhere on your filesystem. For instance in::

  % mkdir -p ~/data/mlcomp
  % cd  ~/data/mlcomp
  % unzip /path/to/dataset-379-20news-18828_XXXXX.zip

You should get a folder ``~/data/mlcomp/379`` with a file named ``metadata`` and
subfolders ``raw``, ``train`` and ``test`` holding the text documents organized 
by newsgroups.

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

from scikits.learn.datasets import load_mlcomp
from scikits.learn.svm.sparse import LinearSVC
from scikits.learn.metrics import confusion_matrix
from scikits.learn.metrics import classification_report

if 'MLCOMP_DATASETS_HOME' not in os.environ:
    print "Please follow those instructions to get started:"
    print __doc__
    sys.exit(0)

# Load the training set
print "Loading 20 newsgroups training set... "
t0 = time()
news_train = load_mlcomp('20news-18828', 'train', sparse=True)
print "done in %fs" % (time() - t0)

print "news_train.data is sparse: ",
print sp.issparse(news_train.data)

# The documents have been hashed into TF-IDF (Term Frequencies times Inverse
# Document Frequencies) vectors of a fixed dimension.
print "n_samples: %d, n_features: %d" % news_train.data.shape

print "Training a linear SVM (hinge loss and L2 regularizer)..."
parameters = {
    'loss': 'l2',
    'penalty': 'l2',
    'C': 10,
    'dual': False,
    'eps': 1e-4,
}
print "parameters:", parameters
t0 = time()
clf = LinearSVC(**parameters).fit(news_train.data, news_train.target)
print "done in %fs" % (time() - t0)
print "Percentage of non zeros coef: %f" % (np.mean(clf.coef_ != 0) * 100)

print "Loading 20 newsgroups test set... "
t0 = time()
news_test = load_mlcomp('20news-18828', 'test', sparse=True)
print "done in %fs" % (time() - t0)

print "Predicting the labels of the test set..."
t0 = time()
pred = clf.predict(news_test.data)
print "done in %fs" % (time() - t0)

print "Classification report on test set:"
print classification_report(news_test.target, pred,
                            class_names=news_test.target_names)


cm = confusion_matrix(news_test.target, pred)
print "Confusion matrix:"
print cm

# Show confusion matrix
pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()
pl.show()
