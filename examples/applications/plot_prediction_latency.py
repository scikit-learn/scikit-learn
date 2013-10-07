"""
==================
Prediction Latency
==================

This is an example showing the prediction latency of various scikit-learn
estimators.
"""

# Authors: Eustache Diemert <eustache@diemert.fr>
# License: BSD 3 clause

from __future__ import print_function

import time
import gc
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import scoreatpercentile
from sklearn.datasets.samples_generator import make_regression
from sklearn.ensemble.forest import RandomForestRegressor
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn.linear_model.ridge import Ridge


def _not_in_sphinx():
    # Hack to detect whether we are running by the sphinx builder
    return '__file__' in globals()


def atomic_benchmark_estimator(estimator, X_test, verbose=False):
    n_instances = X_test.shape[0]
    runtimes = np.zeros(n_instances, dtype=np.float)
    for i in range(n_instances):
        instance = X_test[i, :]
        start = time.time()
        estimator.predict(instance)
        runtimes[i] = time.time() - start
    if verbose:
        print("atomic_benchmark runtimes:", min(runtimes), scoreatpercentile(
            runtimes, 50), max(runtimes))
    return runtimes


def bulk_benchmark_estimator(estimator, X_test, n_bulk_repeats, verbose):
    n_instances = X_test.shape[0]
    runtimes = np.zeros(n_bulk_repeats, dtype=np.float)
    for i in range(n_bulk_repeats):
        start = time.time()
        estimator.predict(X_test)
        runtimes[i] = time.time() - start
    runtimes = np.array(map(lambda x: x/float(n_instances), runtimes))
    if verbose:
        print("bulk_benchmark runtimes:", min(runtimes), scoreatpercentile(
            runtimes, 50), max(runtimes))
    return runtimes


def benchmark_estimator(estimator, X_test, n_bulk_repeats=30, verbose=False):
    """

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.
    """
    atomic_runtimes = atomic_benchmark_estimator(estimator, X_test, verbose)
    bulk_runtimes = bulk_benchmark_estimator(estimator, X_test, n_bulk_repeats,
                                             verbose)
    return atomic_runtimes, bulk_runtimes

def generate_dataset(n_train, n_test, n_features, noise=0.1, verbose=False):
    """
    ...
    """
    if verbose:
        print("generating dataset...")
    X, y, coef = make_regression(n_samples=n_train + n_test,
                                 n_features=n_features, noise=noise, coef=True)
    X_train = X[:n_train]
    y_train = y[:n_train]
    X_test = X[n_train:]
    y_test = y[n_train:]
    idx = np.arange(n_train)
    np.random.seed(13)
    np.random.shuffle(idx)
    X_train = X_train[idx]
    y_train = y_train[idx]

    std = X_train.std(axis=0)
    mean = X_train.mean(axis=0)
    X_train = (X_train - mean) / std
    X_test = (X_test - mean) / std

    std = y_train.std(axis=0)
    mean = y_train.mean(axis=0)
    y_train = (y_train - mean) / std
    y_test = (y_test - mean) / std

    gc.collect()
    if verbose:
        print("ok")
    return X_train, y_train, X_test, y_test

def boxplot_runtimes(runtimes, cls_names, pred_type):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    bp = plt.boxplot(runtimes, )

    xtick_names = plt.setp(ax1, xticklabels=cls_names)
    plt.setp(xtick_names)

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.5)

    ax1.set_axisbelow(True)
    ax1.set_title('Prediction Time per Instance - %s' % pred_type.capitalize())
    ax1.set_xlabel('Estimator')
    ax1.set_ylabel('Prediction Time (us)')

    plt.show()


def benchmark(n_train, n_test, n_feats):
    X_train, y_train, X_test, y_test = generate_dataset(n_train, n_test, n_feats,
                                                        verbose=True)

    stats = {'settings': {'n_train': n_train, 'n_test': 'n_test',
                          'n_feats': n_feats}}
    estimators = {'elasticnet': ElasticNet(), 'ridge': Ridge(),
                  'randomforest': RandomForestRegressor()}
    for clf_name, clf in estimators.iteritems():
        print("Benchmarking", clf)
        clf.fit(X_train, y_train)
        gc.collect()
        a, b = benchmark_estimator(clf, X_test)
        stats[clf_name] = {'atomic': a, 'bulk': b}

    cls_names = estimators.keys()
    runtimes = [1e6*stats[clf_name]['atomic'] for clf_name in cls_names]
    boxplot_runtimes(runtimes, cls_names, ('atomic'))
    runtimes = [1e6*stats[clf_name]['bulk'] for clf_name in cls_names]
    boxplot_runtimes(runtimes, cls_names, 'bulk (%d)' % n_test)

if __name__ == '__main__':
    n_train = int(1e3)
    n_test = int(1e2)
    n_feats = int(1e2)
    benchmark(n_train, n_test, n_feats)
