"""
================================
Classification of text documents
================================

This is an example showing how the scikit-learn can be used to classify
documents by topics using a bag-of-words approach.

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

  % ipython examples/mlcomp_document_classification.py

"""
print __doc__

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

from time import time
import sys
import os
import numpy as np
import pylab as pl

from scikits.learn.datasets import load_mlcomp
from scikits.learn.svm import LinearSVC
from scikits.learn.metrics import confusion_matrix
from scikits.learn.metrics import classification_report

if 'MLCOMP_DATASETS_HOME' not in os.environ:
    print "Please follow those instructions to get started:"
    print __doc__
    sys.exit(0)

# Load the training set
print "Loading 20 newsgroups training set... "
t0 = time()
news_train = load_mlcomp('20news-18828', 'train')
print "done in %fs" % (time() - t0)

# The documents have been hashed into TF-IDF (Term Frequencies times Inverse
# Document Frequencies) vectors of a fixed dimension.
# Currently most scikits.learn wrappers or algorithm implementations are unable
# to leverage efficiently a sparse datastracture; hence we use a dense
# representation of a text dataset. Efficient handling of sparse data
# structures should be expected in an upcoming version of scikits.learn
print "n_samples: %d, n_features: %d" % news_train.data.shape

print "Training a linear classification model with L1 penalty... "
parameters = {
    'loss': 'l1',
    'penalty': 'l2',
    'C': 10,
    'dual': True,
    'eps': 1e-4,
}
print "parameters:", parameters
t0 = time()
clf = LinearSVC(**parameters).fit(news_train.data, news_train.target)
print "done in %fs" % (time() - t0)
print "Percentage of non zeros coef: %f" % (np.mean(clf.coef_ != 0) * 100)

print "Loading 20 newsgroups test set... "
t0 = time()
news_test = load_mlcomp('20news-18828', 'test')
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
