"""
=======================
MNIST dataset benchmark
=======================

Benchmark multi-layer perceptron, Extra-Trees, linear svm
with kernel approximation of RBFSampler and Nystroem
on the MNIST dataset.  The dataset comprises 70,000 samples
and 784 features. Here, we consider the task of predicting
10 classes -  digits from 0 to 9. The experiment was run in
a computer with a Desktop Intel Core i7, 3.6 GHZ CPU,
operating the Windows 7 64-bit version.

    Classification performance:
    ===========================
    Classifier            train-time  test-time error-rate
    ------------------------------------------------------
    MultilayerPerceptron  655.5s      0.30s     0.0169
    nystroem_approx_svm   125.0s      0.91s     0.0239
    ExtraTrees            79.9s       0.34s     0.0272
    fourier_approx_svm    148.9s      0.60s     0.0488
    LogisticRegression    68.9s       0.14s     0.0799


"""


from __future__ import division, print_function

print(__doc__)

# Author: Issam H. Laradji
# License: BSD 3 clause

import logging
import sys
from time import time

import numpy as np
from sklearn import svm, pipeline
from sklearn.datasets import fetch_mldata
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.kernel_approximation import (RBFSampler,
                                          Nystroem)
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MultilayerPerceptronClassifier
from sklearn import metrics


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def load_data(dtype=np.float32, order='F'):
    # Load dataset
    print("Loading dataset...")
    data = fetch_mldata('MNIST original')
    X, y = data.data, data.target

    # Normalize features
    X = X.astype('float64')
    X = X / 255

    # Create train-test split (as [Joachims, 2006])
    print("Creating train-test split...")
    n_train = 60000

    X_train = X[:n_train]
    y_train = y[:n_train]
    X_test = X[n_train:]
    y_test = y[n_train:]

    return X_train, X_test, y_train, y_test


X_train, X_test, y_train, y_test = load_data()

# Print dataset statistics
print("")
print("Dataset statistics:")
print("===================")
print("{0: <25} {1}".format("number of features:",
                            X_train.shape[1]))
print("{0: <25} {1}".format("number of classes:",
                            np.unique(y_train).shape[0]))
print("{0: <25} {1}".format("data type:",
                            X_train.dtype))
print("{0: <25} {1} (pos={2}, neg={3}, size={4}MB)".format(
    "number of train samples:",
      X_train.shape[0], np.sum(y_train == 1),
      np.sum(y_train == -1), int(X_train.nbytes / 1e6)))
print("{0: <25} {1} (pos={2}, neg={3}, size={4}MB)".format(
      "number of test samples:",
      X_test.shape[0], np.sum(y_test == 1),
      np.sum(y_test == -1), int(X_test.nbytes / 1e6)))


classifiers = dict()


# Benchmark classifiers
def benchmark(clf):
    t0 = time()
    clf.fit(X_train, y_train)
    train_time = time() - t0
    t0 = time()
    pred = clf.predict(X_test)
    test_time = time() - t0
    err = metrics.zero_one_loss(y_test, pred)
    return err, train_time, test_time

# Train MultilayerPerceptron model
classifiers['LogisticRegression'] = LogisticRegression()

# Train MultilayerPerceptron model
classifiers['MultilayerPerceptron'] = MultilayerPerceptronClassifier(
    hidden_layer_sizes=(200,),
    max_iter=400,
    alpha=0.5,
    algorithm='l-bfgs',
    random_state=1)

# Train Extra-Trees model
classifiers['ExtraTrees'] = ExtraTreesClassifier(n_estimators=100,
                                                 random_state=1)

# Train linear svm with kernel approximation of RBFSampler and Nystroem

# create pipeline
feature_map_fourier = RBFSampler(gamma=.031, random_state=1,
                                 n_components=1000)
feature_map_nystroem = Nystroem(gamma=.031, random_state=1, n_components=1000)

classifiers['fourier_approx_svm'] = \
    pipeline.Pipeline(
        [("feature_map",
          feature_map_fourier),
         ("svm", svm.LinearSVC(C=100))])

classifiers['nystroem_approx_svm'] = \
    pipeline.Pipeline([("feature_map",
                        feature_map_nystroem),
                       ("svm", svm.LinearSVC(C=100, random_state=1))])

selected_classifiers = classifiers.keys()
for name in selected_classifiers:
    if name not in classifiers:
        print('classifier %r unknown' % name)
        sys.exit(1)

print()
print("Training Classifiers")
print("====================")
print()
err, train_time, test_time = {}, {}, {}
for name in sorted(selected_classifiers):
    print("Training {0} ...".format(name))
    err[name], train_time[name], test_time[name] = benchmark(classifiers[name])

# Print classification performance
print()
print("Classification performance:")
print("===========================")
print()


def print_row(clf_type, train_time, test_time, err):
    print("{0: <25} {1: ^15} {2: ^15} {3: ^15}".format(
          clf_type,
          train_time,
          test_time,
          err))

print("{0: <25} {1: ^15} {2: ^15} {3: ^15}".format(
    "Classifier  ",
    "train-time",
    "test-time",
    "error-rate"))
print("-" * 70)

for name in sorted(selected_classifiers, key=lambda name: err[name]):
    print_row(name, train_time[name], test_time[name], err[name])
print()
print()
