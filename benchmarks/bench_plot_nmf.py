"""
Benchmarks of Non-Negative Matrix Factorization
"""
# Author : Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
# License: BSD 3 clause

from __future__ import print_function
from time import time
import sys

import six

import numpy as np
import matplotlib.pyplot as plt
import pandas

from sklearn.utils.testing import ignore_warnings
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition.nmf import NMF
from sklearn.decomposition.nmf import _initialize_nmf
from sklearn.decomposition.nmf import _beta_divergence
from sklearn.externals.joblib import Memory
from sklearn.exceptions import ConvergenceWarning

mem = Memory(cachedir='.', verbose=0)


def plot_results(results_df, plot_name):
    if results_df is None:
        return None

    plt.figure(figsize=(16, 6))
    colors = 'bgr'
    markers = 'ovs'
    ax = plt.subplot(1, 3, 1)
    for i, init in enumerate(np.unique(results_df['init'])):
        plt.subplot(1, 3, i + 1, sharex=ax, sharey=ax)
        for j, method in enumerate(np.unique(results_df['method'])):
            mask = np.logical_and(results_df['init'] == init,
                                  results_df['method'] == method)
            selected_items = results_df[mask]

            plt.plot(selected_items['time'], selected_items['loss'],
                     color=colors[j % len(colors)], ls='-',
                     marker=markers[j % len(markers)],
                     label=method)

        plt.legend(loc=0, fontsize='x-small')
        plt.xlabel("Time (s)")
        plt.ylabel("loss")
        plt.title("%s" % init)
    plt.suptitle(plot_name, fontsize=16)


# The deprecated projected-gradient solver raises a UserWarning as convergence
# is not reached; the coordinate-descent solver raises a ConvergenceWarning.
@ignore_warnings(category=(ConvergenceWarning, UserWarning,
                           DeprecationWarning))
# use joblib to cache the results.
# X_shape is specified in arguments for avoiding hashing X
@mem.cache(ignore=['X', 'W0', 'H0'])
def bench_one(name, X, W0, H0, X_shape, clf_type, clf_params, init,
              n_components, random_state):
    W = W0.copy()
    H = H0.copy()

    clf = clf_type(**clf_params)
    st = time()
    W = clf.fit_transform(X, W=W, H=H)
    end = time()
    H = clf.components_

    this_loss = _beta_divergence(X, W, H, 2.0, True)
    duration = end - st
    return this_loss, duration


def run_bench(X, clfs, plot_name, n_components, tol, alpha, l1_ratio):
    start = time()
    results = []
    for name, clf_type, iter_range, clf_params in clfs:
        print("Training %s:" % name)
        for rs, init in enumerate(('nndsvd', 'nndsvdar', 'random')):
            print("    %s %s: " % (init, " " * (8 - len(init))), end="")
            W, H = _initialize_nmf(X, n_components, init, 1e-6, rs)

            for max_iter in iter_range:
                clf_params['alpha'] = alpha
                clf_params['l1_ratio'] = l1_ratio
                clf_params['max_iter'] = max_iter
                clf_params['tol'] = tol
                clf_params['random_state'] = rs
                clf_params['init'] = 'custom'
                clf_params['n_components'] = n_components

                this_loss, duration = bench_one(name, X, W, H, X.shape,
                                                clf_type, clf_params,
                                                init, n_components, rs)

                init_name = "init='%s'" % init
                results.append((name, this_loss, duration, init_name))
                # print("loss: %.6f, time: %.3f sec" % (this_loss, duration))
                print(".", end="")
                sys.stdout.flush()
            print(" ")

    # Use a panda dataframe to organize the results
    results_df = pandas.DataFrame(results,
                                  columns="method loss time init".split())
    print("Total time = %0.3f sec\n" % (time() - start))

    # plot the results
    plot_results(results_df, plot_name)
    return results_df


def load_20news():
    print("Loading 20 newsgroups dataset")
    print("-----------------------------")
    from sklearn.datasets import fetch_20newsgroups
    dataset = fetch_20newsgroups(shuffle=True, random_state=1,
                                 remove=('headers', 'footers', 'quotes'))
    vectorizer = TfidfVectorizer(max_df=0.95, min_df=2, stop_words='english')
    tfidf = vectorizer.fit_transform(dataset.data)
    return tfidf


def load_faces():
    print("Loading Olivetti face dataset")
    print("-----------------------------")
    from sklearn.datasets import fetch_olivetti_faces
    faces = fetch_olivetti_faces(shuffle=True)
    return faces.data


def build_clfs(cd_iters, mu_iters):
    clfs = [("Coordinate Descent", NMF, cd_iters, {'solver': 'cd'}),
            ("Multiplicative Update", NMF, mu_iters, {'solver': 'mu'}),
            ]
    return clfs


if __name__ == '__main__':
    alpha = 0.
    l1_ratio = 0.5
    n_components = 10
    tol = 1e-15

    # first benchmark on 20 newsgroup dataset: sparse, shape(11314, 39116)
    plot_name = "20 Newsgroups sparse dataset"
    cd_iters = np.arange(1, 30)
    mu_iters = np.arange(1, 30)
    clfs = build_clfs(cd_iters, mu_iters)
    X_20news = load_20news()
    run_bench(X_20news, clfs, plot_name, n_components, tol, alpha, l1_ratio)

    # second benchmark on Olivetti faces dataset: dense, shape(400, 4096)
    plot_name = "Olivetti Faces dense dataset"
    cd_iters = np.arange(1, 30)
    mu_iters = np.arange(1, 30)
    clfs = build_clfs(cd_iters, mu_iters)
    X_faces = load_faces()
    run_bench(X_faces, clfs, plot_name, n_components, tol, alpha, l1_ratio,)

    plt.show()
