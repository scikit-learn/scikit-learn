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


def atomic_benchmark_estimator(estimator, X_test, percentile, verbose=False):
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
    return scoreatpercentile(runtimes, percentile)


def bulk_benchmark_estimator(estimator, X_test, percentile, n_bulk_repeats,
                             verbose):
    n_instances = X_test.shape[0]
    runtimes = np.zeros(n_bulk_repeats, dtype=np.float)
    for i in range(n_bulk_repeats):
        start = time.time()
        estimator.predict(X_test)
        runtimes[i] = time.time() - start
    runtimes = map(lambda x: x/float(n_instances), runtimes)
    if verbose:
        print("bulk_benchmark runtimes:", min(runtimes), scoreatpercentile(
            runtimes, 50), max(runtimes))
    return scoreatpercentile(runtimes, percentile)


def benchmark_estimator(estimator, X_test, percentile=90, n_bulk_repeats=30,
                        verbose=False):
    """

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.
    """
    atomic_latency = atomic_benchmark_estimator(estimator, X_test,
                                                percentile, verbose)
    bulk_latency = bulk_benchmark_estimator(estimator, X_test,
                                            percentile, n_bulk_repeats,
                                            verbose)
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
n_feats = int(1e2)
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

plt.figure()
ax = plt.subplot(111)
cls_names = estimators.keys()
runtimes = [1e6*stats[clf_name]['atomic'] for clf_name in cls_names]
bar_colors = rcParams['axes.color_cycle'][:len(cls_names)]
rectangles = plt.bar(range(len(stats)-1),
                     runtimes,
                     width=0.5,
                     color=bar_colors)

ax.set_xticks(np.linspace(0.25, len(cls_names) - 0.75, len(cls_names)))
ax.set_xticklabels(cls_names, fontsize=10)
ymax = max(runtimes) * 1.2
ax.set_ylim((0, ymax))
ax.set_ylabel('runtime (micro-seconds) per instance')
ax.set_title('Prediction Latency (atomic)')


def autolabel(rectangles):
    """attach some text vi autolabel on rectangles."""
    for rect in rectangles:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2.,
                1.05 * height, '%.4f' % height,
                ha='center', va='bottom')

autolabel(rectangles)
plt.show()
