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
   <https://www.cs.cornell.edu/people/tj/publications/joachims_06a.pdf>`_
   T. Joachims - In SIGKDD '06

[1] https://archive.ics.uci.edu/ml/datasets/Covertype

"""
from __future__ import division, print_function

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Arnaud Joly <arnaud.v.joly@gmail.com>
# License: BSD 3 clause

import os
from time import time
import argparse
import numpy as np
from joblib import Memory

from sklearn.datasets import fetch_covtype, get_data_home
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import zero_one_loss
from sklearn.utils import check_array

# Memoize the data extraction and memory map the resulting
# train / test splits in readonly mode
memory = Memory(os.path.join(get_data_home(), 'covertype_benchmark_data'),
                mmap_mode='r')


@memory.cache
def load_data(dtype=np.float32, order='C', random_state=13):
    """Load the data, then cache and memmap the train/test split"""
    ######################################################################
    # Load dataset
    print("Loading dataset...")
    data = fetch_covtype(download_if_missing=True, shuffle=True,
                         random_state=random_state)
    X = check_array(data['data'], dtype=dtype, order=order)
    y = (data['target'] != 1).astype(np.int)

    # Create train-test split (as [Joachims, 2006])
    print("Creating train-test split...")
    n_train = 522911
    X_train = X[:n_train]
    y_train = y[:n_train]
    X_test = X[n_train:]
    y_test = y[n_train:]

    # Standardize first 10 features (the numerical ones)
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    mean[10:] = 0.0
    std[10:] = 1.0
    X_train = (X_train - mean) / std
    X_test = (X_test - mean) / std
    return X_train, X_test, y_train, y_test


ESTIMATORS = {
    'GBRT': GradientBoostingClassifier(n_estimators=250),
    'ExtraTrees': ExtraTreesClassifier(n_estimators=20),
    'RandomForest': RandomForestClassifier(n_estimators=20),
    'CART': DecisionTreeClassifier(min_samples_split=5),
    'SGD': SGDClassifier(alpha=0.001, max_iter=1000, tol=1e-3),
    'GaussianNB': GaussianNB(),
    'liblinear': LinearSVC(loss="l2", penalty="l2", C=1000, dual=False,
                           tol=1e-3),
    'SAG': LogisticRegression(solver='sag', max_iter=2, C=1000)
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--classifiers', nargs="+",
                        choices=ESTIMATORS, type=str,
                        default=['liblinear', 'GaussianNB', 'SGD', 'CART'],
                        help="list of classifiers to benchmark.")
    parser.add_argument('--n-jobs', nargs="?", default=1, type=int,
                        help="Number of concurrently running workers for "
                             "models that support parallelism.")
    parser.add_argument('--order', nargs="?", default="C", type=str,
                        choices=["F", "C"],
                        help="Allow to choose between fortran and C ordered "
                             "data")
    parser.add_argument('--random-seed', nargs="?", default=13, type=int,
                        help="Common seed used by random number generator.")
    args = vars(parser.parse_args())

    print(__doc__)

    X_train, X_test, y_train, y_test = load_data(
        order=args["order"], random_state=args["random_seed"])

    print("")
    print("Dataset statistics:")
    print("===================")
    print("%s %d" % ("number of features:".ljust(25), X_train.shape[1]))
    print("%s %d" % ("number of classes:".ljust(25), np.unique(y_train).size))
    print("%s %s" % ("data type:".ljust(25), X_train.dtype))
    print("%s %d (pos=%d, neg=%d, size=%dMB)"
          % ("number of train samples:".ljust(25),
             X_train.shape[0], np.sum(y_train == 1),
             np.sum(y_train == 0), int(X_train.nbytes / 1e6)))
    print("%s %d (pos=%d, neg=%d, size=%dMB)"
          % ("number of test samples:".ljust(25),
             X_test.shape[0], np.sum(y_test == 1),
             np.sum(y_test == 0), int(X_test.nbytes / 1e6)))

    print()
    print("Training Classifiers")
    print("====================")
    error, train_time, test_time = {}, {}, {}
    for name in sorted(args["classifiers"]):
        print("Training %s ... " % name, end="")
        estimator = ESTIMATORS[name]
        estimator_params = estimator.get_params()

        estimator.set_params(**{p: args["random_seed"]
                                for p in estimator_params
                                if p.endswith("random_state")})

        if "n_jobs" in estimator_params:
            estimator.set_params(n_jobs=args["n_jobs"])

        time_start = time()
        estimator.fit(X_train, y_train)
        train_time[name] = time() - time_start

        time_start = time()
        y_pred = estimator.predict(X_test)
        test_time[name] = time() - time_start

        error[name] = zero_one_loss(y_test, y_pred)

        print("done")

    print()
    print("Classification performance:")
    print("===========================")
    print("%s %s %s %s"
          % ("Classifier  ", "train-time", "test-time", "error-rate"))
    print("-" * 44)
    for name in sorted(args["classifiers"], key=error.get):
        print("%s %s %s %s" % (name.ljust(12),
                               ("%.4fs" % train_time[name]).center(10),
                               ("%.4fs" % test_time[name]).center(10),
                               ("%.4f" % error[name]).center(10)))

    print()
