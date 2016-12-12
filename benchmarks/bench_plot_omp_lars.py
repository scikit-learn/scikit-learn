"""Benchmarks of orthogonal matching pursuit (:ref:`OMP`) versus least angle
regression (:ref:`least_angle_regression`)

The input data is mostly low rank but is a fat infinite tail.
"""
from __future__ import print_function

import gc
import sys
from time import time

import six

import numpy as np

from sklearn.linear_model import lars_path, orthogonal_mp
from sklearn.datasets.samples_generator import make_sparse_coded_signal


def compute_bench(samples_range, features_range):

    it = 0

    results = dict()
    lars = np.empty((len(features_range), len(samples_range)))
    lars_gram = lars.copy()
    omp = lars.copy()
    omp_gram = lars.copy()

    max_it = len(samples_range) * len(features_range)
    for i_s, n_samples in enumerate(samples_range):
        for i_f, n_features in enumerate(features_range):
            it += 1
            n_informative = n_features / 10
            print('====================')
            print('Iteration %03d of %03d' % (it, max_it))
            print('====================')
            # dataset_kwargs = {
            #     'n_train_samples': n_samples,
            #     'n_test_samples': 2,
            #     'n_features': n_features,
            #     'n_informative': n_informative,
            #     'effective_rank': min(n_samples, n_features) / 10,
            #     #'effective_rank': None,
            #     'bias': 0.0,
            # }
            dataset_kwargs = {
                'n_samples': 1,
                'n_components': n_features,
                'n_features': n_samples,
                'n_nonzero_coefs': n_informative,
                'random_state': 0
            }
            print("n_samples: %d" % n_samples)
            print("n_features: %d" % n_features)
            y, X, _ = make_sparse_coded_signal(**dataset_kwargs)
            X = np.asfortranarray(X)

            gc.collect()
            print("benchmarking lars_path (with Gram):", end='')
            sys.stdout.flush()
            tstart = time()
            G = np.dot(X.T, X)  # precomputed Gram matrix
            Xy = np.dot(X.T, y)
            lars_path(X, y, Xy=Xy, Gram=G, max_iter=n_informative)
            delta = time() - tstart
            print("%0.3fs" % delta)
            lars_gram[i_f, i_s] = delta

            gc.collect()
            print("benchmarking lars_path (without Gram):", end='')
            sys.stdout.flush()
            tstart = time()
            lars_path(X, y, Gram=None, max_iter=n_informative)
            delta = time() - tstart
            print("%0.3fs" % delta)
            lars[i_f, i_s] = delta

            gc.collect()
            print("benchmarking orthogonal_mp (with Gram):", end='')
            sys.stdout.flush()
            tstart = time()
            orthogonal_mp(X, y, precompute=True,
                          n_nonzero_coefs=n_informative)
            delta = time() - tstart
            print("%0.3fs" % delta)
            omp_gram[i_f, i_s] = delta

            gc.collect()
            print("benchmarking orthogonal_mp (without Gram):", end='')
            sys.stdout.flush()
            tstart = time()
            orthogonal_mp(X, y, precompute=False,
                          n_nonzero_coefs=n_informative)
            delta = time() - tstart
            print("%0.3fs" % delta)
            omp[i_f, i_s] = delta

    results['time(LARS) / time(OMP)\n (w/ Gram)'] = (lars_gram / omp_gram)
    results['time(LARS) / time(OMP)\n (w/o Gram)'] = (lars / omp)
    return results


if __name__ == '__main__':
    samples_range = np.linspace(1000, 5000, 5).astype(np.int)
    features_range = np.linspace(1000, 5000, 5).astype(np.int)
    results = compute_bench(samples_range, features_range)
    max_time = max(np.max(t) for t in results.values())

    import matplotlib.pyplot as plt
    fig = plt.figure('scikit-learn OMP vs. LARS benchmark results')
    for i, (label, timings) in enumerate(sorted(six.iteritems(results))):
        ax = fig.add_subplot(1, 2, i+1)
        vmax = max(1 - timings.min(), -1 + timings.max())
        plt.matshow(timings, fignum=False, vmin=1 - vmax, vmax=1 + vmax)
        ax.set_xticklabels([''] + [str(each) for each in samples_range])
        ax.set_yticklabels([''] + [str(each) for each in features_range])
        plt.xlabel('n_samples')
        plt.ylabel('n_features')
        plt.title(label)

    plt.subplots_adjust(0.1, 0.08, 0.96, 0.98, 0.4, 0.63)
    ax = plt.axes([0.1, 0.08, 0.8, 0.06])
    plt.colorbar(cax=ax, orientation='horizontal')
    plt.show()
