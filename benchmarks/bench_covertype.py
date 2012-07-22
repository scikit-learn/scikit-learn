"""
===========================
Covertype dataset benchmark
===========================

Benchmark stochastic gradient descent (SGD), Liblinear, and Naive Bayes, CART
(decision tree), RandomForest and Extra-Trees on the forest covertype dataset
of Blackard, Jock, and Dean [1]. The dataset comprises 581,012 samples. It is
low dimensional with 54 features and a sparsity of approx. 23%. Here, we
consider the task of predicting class 1 (spruce/fir). The classification
performance of SGD is competitive with Liblinear while being two orders of
magnitude faster to train::

    [..]
    Classification performance:
    ===========================

    Classifier   train-time test-time error-rate
    --------------------------------------------
    Liblinear     11.8977s   0.0285s     0.2305
    GaussianNB    3.5931s    0.6645s     0.3633
    SGD           0.2924s    0.0114s     0.2300
    CART          39.9829s   0.0345s     0.0476
    RandomForest  794.6232s  1.0526s     0.0249
    Extra-Trees   1401.7051s 1.1181s     0.0230

The same task has been used in a number of papers including:

 * `"SVM Optimization: Inverse Dependence on Training Set Size"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.139.2112>`_
   S. Shalev-Shwartz, N. Srebro - In Proceedings of ICML '08.

 * `"Pegasos: Primal estimated sub-gradient solver for svm"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.74.8513>`_
   S. Shalev-Shwartz, Y. Singer, N. Srebro - In Proceedings of ICML '07.

 * `"Training Linear SVMs in Linear Time"
   <www.cs.cornell.edu/People/tj/publications/joachims_06a.pdf>`_
   T. Joachims - In SIGKDD '06

[1] http://archive.ics.uci.edu/ml/datasets/Covertype

"""
from __future__ import division

print __doc__

# Author: Peter Prettenhoer <peter.prettenhofer@gmail.com>
# License: BSD Style.

# $Id$

from time import time
import os
import sys
import numpy as np
from optparse import OptionParser

from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn import metrics

op = OptionParser()
op.add_option("--classifiers",
              dest="classifiers", default='liblinear,GaussianNB,SGD,CART',
              help="comma-separated list of classifiers to benchmark. "
                   "default: %default. available: "
                   "liblinear,GaussianNB,SGD,CART,ExtraTrees,RandomForest")

op.print_help()

(opts, args) = op.parse_args()
if len(args) > 0:
    op.error("this script takes no arguments.")
    sys.exit(1)

######################################################################
## Download the data, if not already on disk
if not os.path.exists('covtype.data.gz'):
    # Download the data
    import urllib
    print "Downloading data, Please Wait (11MB)..."
    opener = urllib.urlopen(
        'http://archive.ics.uci.edu/ml/'
        'machine-learning-databases/covtype/covtype.data.gz')
    open('covtype.data.gz', 'wb').write(opener.read())

######################################################################
## Load dataset
print("Loading dataset...")
import gzip
f = gzip.open('covtype.data.gz')
X = np.fromstring(f.read().replace(",", " "), dtype=np.float64, sep=" ",
                  count=-1)
X = X.reshape((581012, 55))
f.close()

# class 1 vs. all others.
y = np.ones(X.shape[0]) * -1
y[np.where(X[:, -1] == 1)] = 1
X = X[:, :-1]

######################################################################
## Create train-test split (as [Joachims, 2006])
print("Creating train-test split...")
idx = np.arange(X.shape[0])
np.random.seed(13)
np.random.shuffle(idx)
train_idx = idx[:522911]
test_idx = idx[522911:]

X_train = X[train_idx]
y_train = y[train_idx]
X_test = X[test_idx]
y_test = y[test_idx]

# free memory
del X
del y

######################################################################
## Standardize first 10 features (the numerical ones)
mean = X_train.mean(axis=0)
std = X_train.std(axis=0)
mean[10:] = 0.0
std[10:] = 1.0
X_train = (X_train - mean) / std
X_test = (X_test - mean) / std

######################################################################
## Print dataset statistics
print("")
print("Dataset statistics:")
print("===================")
print("%s %d" % ("number of features:".ljust(25),
                 X_train.shape[1]))
print("%s %d" % ("number of classes:".ljust(25),
                 np.unique(y_train).shape[0]))
print("%s %d (%d, %d)" % ("number of train samples:".ljust(25),
                          X_train.shape[0], np.sum(y_train == 1),
                          np.sum(y_train == -1)))
print("%s %d (%d, %d)" % ("number of test samples:".ljust(25),
                          X_test.shape[0], np.sum(y_test == 1),
                          np.sum(y_test == -1)))


classifiers = dict()


######################################################################
## Benchmark classifiers
def benchmark(clf):
    t0 = time()
    clf.fit(X_train, y_train)
    train_time = time() - t0
    t0 = time()
    pred = clf.predict(X_test)
    test_time = time() - t0
    err = metrics.zero_one(y_test, pred) / float(pred.shape[0])
    return err, train_time, test_time

######################################################################
## Train Liblinear model
liblinear_parameters = {
    'loss': 'l2',
    'penalty': 'l2',
    'C': 1000,
    'dual': False,
    'tol': 1e-3,
    }
classifiers['liblinear'] = LinearSVC(**liblinear_parameters)

######################################################################
## Train GaussianNB model
classifiers['GaussianNB'] = GaussianNB()

######################################################################
## Train SGD model
sgd_parameters = {
    'alpha': 0.001,
    'n_iter': 2,
    }
classifiers['SGD'] = SGDClassifier(**sgd_parameters)

######################################################################
## Train CART model
classifiers['CART'] = DecisionTreeClassifier(min_samples_split=5,
                                             max_depth=None)

######################################################################
## Train RandomForest model
classifiers['RandomForest'] = RandomForestClassifier(n_estimators=20,
                                                     min_samples_split=5,
                                                     max_features=None,
                                                     max_depth=None)

######################################################################
## Train Extra-Trees model
classifiers['ExtraTrees'] = ExtraTreesClassifier(n_estimators=20,
                                                 min_samples_split=5,
                                                 max_features=None,
                                                 max_depth=None)


selected_classifiers = opts.classifiers.split(',')
for name in selected_classifiers:
    if name not in classifiers:
        op.error('classifier %r unknown' % name)
        sys.exit(1)

print("")
print("Training Classifiers")
print("====================")
print("")
err, train_time, test_time = {}, {}, {}
for name in sorted(selected_classifiers):
    print("Training %s ..." % name)
    err[name], train_time[name], test_time[name] = benchmark(classifiers[name])

######################################################################
## Print classification performance
print("")
print("Classification performance:")
print("===========================")
print("")


def print_row(clf_type, train_time, test_time, err):
    print("%s %s %s %s" % (clf_type.ljust(12),
                           ("%.4fs" % train_time).center(10),
                           ("%.4fs" % test_time).center(10),
                           ("%.4f" % err).center(10)))

print("%s %s %s %s" % ("Classifier  ", "train-time", "test-time",
                       "error-rate"))
print("-" * 44)

for name in sorted(selected_classifiers, key=lambda name: err[name]):
    print_row(name, train_time[name], test_time[name], err[name])
print("")
print("")
