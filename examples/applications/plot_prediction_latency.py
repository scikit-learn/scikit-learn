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
import numpy as np

from scipy.stats import scoreatpercentile
from sklearn.datasets.samples_generator import make_regression
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn.linear_model.ridge import Ridge


def _not_in_sphinx():
    # Hack to detect whether we are running by the sphinx builder
    return '__file__' in globals()

def benchmark_estimator(estimator, X_test, percentile=90, bulk_repeats=30,
                        verbose=False):
    """

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.
    """
    runtimes = []
    for instanceidx in range(X_test.shape[0]):
        start = time.time()
        estimator.predict(X_test[instanceidx, :])
        runtimes.append(time.time() - start)
    atomic_latency = scoreatpercentile(np.array(runtimes), percentile)
    runtimes = []
    for _repeat in range(bulk_repeats):
        start = time.time()
        clf.predict(X_test)
        runtimes.append(time.time() - start)
    runtimes = map(lambda x:x/float(bulk_repeats), runtimes)
    bulk_latency = scoreatpercentile(np.array(runtimes), percentile)
    if verbose:
        print("prediction time (atomic) <= %.3f us/instance at %d%%" % (
            atomic_latency*1e6, percentile))
        print("prediction time (bulk) <= %.3f us/instance at %d%%" % (
            bulk_latency*1e6, percentile))
    return atomic_latency, bulk_latency

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

n_train = int(1e4)
n_test = int(1e3)
n_feats = int(1e3)
X_train, y_train, X_test, y_test = generate_dataset(n_train, n_test, n_feats,
                                                    verbose=True)
for clf in (ElasticNet(alpha=0.01, l1_ratio=0.5, fit_intercept=False),
            Ridge(alpha=0.01, fit_intercept=False)):
    print(clf)
    clf.fit(X_train, y_train)
    benchmark_estimator(clf, X_test, verbose=True)

