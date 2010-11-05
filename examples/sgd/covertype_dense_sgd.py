"""
================================
Covertype dataset with dense SGD
================================

This is an example showing how stochastic gradient descent (SGD)
can be used for large and dense datasets. 

The dataset used in this example is the covertype dataset of Blackard, Jock, and Dean, which is low-dimensional with 54 features and a sparsity of 23%. The dataset
comprises 581,012 samples. The task is to predicting forest cover
type from cartographic variables only. Each sample represents a 30x30 meter cell
and there are 7 types of forest cover. See [1] for further details.
In this example, however, we consider the task of discriminating between class 1
and all others.

The same task has been used in a number of papers including:

 * `"Pegasos: Primal estimated sub-gradient solver for svm" 
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.74.8513>`_
   S. Shalev-Shwartz, Y. Singer, N. Srebro - In Proceedings of ICML '07.

 * "Training Linear SVMs in Linear Time"  
   T. Joachims - In SIGKDD '06


[1] http://archive.ics.uci.edu/ml/datasets/Covertype

To run this example use your favorite python shell::

  % ipython examples/sgd/covertype_dense_sgd.py

"""
from __future__ import division

print __doc__

# Author: Peter Prettenhoer <peter.prettenhofer@gmail.com>
#
# License: Simplified BSD

from time import time
import sys
import os
import numpy as np
from scipy import sparse
import pylab as pl

from scikits.learn.svm import LinearSVC
from scikits.learn.sgd import SGD
from scikits.learn import metrics

# Download the data, if not already on disk 
if not os.path.exists('covtype.data.gz'):
    # Download the data
    import urllib
    print "Downloading data, Please Wait (11MB)"
    opener = urllib.urlopen(
        'http://archive.ics.uci.edu/ml/' \
        'machine-learning-databases/covtype/covtype.data.gz')
    open('covtype.data.gz', 'wb').write(opener.read())

# Load dataset
print("loading dataset...")
data = np.loadtxt('covtype.data.gz', delimiter=",")
X = data[:, :-1]

# class 1 vs. all others. 
y = np.ones(data.shape[0]) * -1.0
y[np.where(data[:, -1] == 1)] = 1.0

print("creating train-test split...")
idx = np.arange(data.shape[0])
np.random.seed(13)
np.random.shuffle(idx)
train_idx = idx[:522911]
test_idx = idx[522911:]

X_train = X[train_idx]
y_train = y[train_idx]

X_test = X[test_idx]
y_test = y[test_idx]

# standardize first 10 features.
mean = X_train.mean(axis=0)
std = X_train.std(axis=0)
mean[10:] = 0.0
std[10:] = 1.0

X_train = (X_train-mean) / std
X_test = (X_test-mean) / std

# If you use sparse.SGD it's even faster
#X_train = sparse.csr_matrix(X_train)
#X_test = sparse.csr_matrix(X_test)

print("num training: %d" % train_idx.shape[0])
print("num test: %d" % test_idx.shape[0])

## parameters = {
##     'loss': 'l2',
##     'penalty': 'l2',
##     'C': 1000,
##     'dual': False,
##     'eps': 1e-3,
## }
## print("Training LinearSVC on training set")
## clf = LinearSVC(**parameters)
print("Training SGD with alpha=0.001 and n_iter=2")
clf = SGD(alpha=0.001, n_iter=2)
t0 = time()
clf.fit(X_train, y_train)
print "done in %fs" % (time() - t0)

print "Predicting the outcomes of the testing set"
t0 = time()
pred = clf.predict(X_test)
print "done in %fs" % (time() - t0)

print "Classification performance:"
print
err = metrics.zero_one(y_test, pred) / float(pred.shape[0])
print "Errorrate: %.4f" % err
cm = metrics.confusion_matrix(y_test, pred)
print "Confusion matrix:"
print cm

# Show confusion matrix
pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()
pl.show()


