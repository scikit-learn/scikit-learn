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
    liblinear     15.9744s    0.0705s     0.2305
    GaussianNB    3.0666s     0.3884s     0.4841
    SGD           1.0558s     0.1152s     0.2300
    CART          79.4296s    0.0523s     0.0469
    RandomForest  1190.1620s  0.5881s     0.0243
    ExtraTrees    640.3194s   0.6495s     0.0198

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
from __future__ import division, print_function

print(__doc__)

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

import logging
import os
import sys
from time import time
from optparse import OptionParser

import numpy as np

from sklearn.datasets import fetch_covtype, get_data_home
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import metrics
from sklearn.externals.joblib import Memory

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

op = OptionParser()
op.add_option("--classifiers",
              dest="classifiers", default='liblinear,GaussianNB,SGD,CART',
              help="comma-separated list of classifiers to benchmark. "
                   "default: %default. available: "
                   "liblinear, GaussianNB, SGD, CART, ExtraTrees,\n"
                   "RandomForest, GBRT")

op.add_option("--n-jobs",
              dest="n_jobs", default=1, type=int,
              help="Number of concurrently running workers for models that"
                   " support parallelism.")

# Each number generator use the same seed to avoid coupling issue between
# estimators.
op.add_option("--random-seed",
              dest="random_seed", default=13, type=int,
              help="Common seed used by random number generator.")

op.print_help()

(opts, args) = op.parse_args()
if len(args) > 0:
    op.error("this script takes no arguments.")
    sys.exit(1)

# Memoize the data extraction and memory map the resulting
# train / test splits in readonly mode
joblib_cache_folder = os.path.join(get_data_home(), 'covertype_benchmark_data')
m = Memory(joblib_cache_folder, mmap_mode='r')


# Load the data, then cache and memmap the train/test split
@m.cache
def load_data(dtype=np.float32, order='C'):
    ######################################################################
    ## Load dataset
    print("Loading dataset...")
    data = fetch_covtype(download_if_missing=True, shuffle=True,
                         random_state=opts.random_seed)
    X, y = data['data'], data['target']
    X = np.asarray(X, dtype=dtype)
    
    if order.lower() == 'f':
        X = np.asfortranarray(X)

    # class 1 vs. all others.
    y[np.where(y != 1)] = -1

    ######################################################################
    ## Create train-test split (as [Joachims, 2006])
    logger.info("Creating train-test split...")
    n_train = 522911

    X_train = X[:n_train]
    y_train = y[:n_train]
    X_test = X[n_train:]
    y_test = y[n_train:]

    ######################################################################
    ## Standardize first 10 features (the numerical ones)
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    mean[10:] = 0.0
    std[10:] = 1.0
    X_train = (X_train - mean) / std
    X_test = (X_test - mean) / std
    return X_train, X_test, y_train, y_test


X_train, X_test, y_train, y_test = load_data()

######################################################################
## Print dataset statistics
print("")
print("Dataset statistics:")
print("===================")
print("%s %d" % ("number of features:".ljust(25),
                 X_train.shape[1]))
print("%s %d" % ("number of classes:".ljust(25),
                 np.unique(y_train).shape[0]))
print("%s %s" % ("data type:".ljust(25), X_train.dtype))
print("%s %d (pos=%d, neg=%d, size=%dMB)"
      % ("number of train samples:".ljust(25),
         X_train.shape[0], np.sum(y_train == 1),
         np.sum(y_train == -1), int(X_train.nbytes / 1e6)))
print("%s %d (pos=%d, neg=%d, size=%dMB)"
      % ("number of test samples:".ljust(25),
      X_test.shape[0], np.sum(y_test == 1),
      np.sum(y_test == -1), int(X_test.nbytes / 1e6)))


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
    err = metrics.zero_one_loss(y_test, pred, normalize=True)
    return err, train_time, test_time

######################################################################
## Train Liblinear model
liblinear_parameters = {
    'loss': 'l2',
    'penalty': 'l2',
    'C': 1000,
    'dual': False,
    'tol': 1e-3,
    "random_state": opts.random_seed,
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
    'n_jobs': opts.n_jobs,
    "random_state": opts.random_seed,
}
classifiers['SGD'] = SGDClassifier(**sgd_parameters)

######################################################################
## Train CART model
classifiers['CART'] = DecisionTreeClassifier(min_samples_split=5,
                                             max_depth=None,
                                             random_state=opts.random_seed)

######################################################################
## Train RandomForest model
rf_parameters = {
    "n_estimators": 20,
    "n_jobs": opts.n_jobs,
    "random_state": opts.random_seed,
}
classifiers['RandomForest'] = RandomForestClassifier(**rf_parameters)

######################################################################
## Train Extra-Trees model
classifiers['ExtraTrees'] = ExtraTreesClassifier(n_estimators=20,
                                                 n_jobs=opts.n_jobs,
                                                 random_state=opts.random_seed)

######################################################################
## Train GBRT model
classifiers['GBRT'] = GradientBoostingClassifier(n_estimators=250,
                                                 random_state=opts.random_seed)


selected_classifiers = opts.classifiers.split(',')
for name in selected_classifiers:
    if name not in classifiers:
        op.error('classifier %r unknown' % name)
        sys.exit(1)

print()
print("Training Classifiers")
print("====================")
print()
err, train_time, test_time = {}, {}, {}
for name in sorted(selected_classifiers):
    print("Training %s ..." % name)
    err[name], train_time[name], test_time[name] = benchmark(classifiers[name])

######################################################################
## Print classification performance
print()
print("Classification performance:")
print("===========================")
print()


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
print()
print()
