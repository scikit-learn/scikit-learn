"""
===========================
Random projection benchmark
===========================

Transformer      |  dense-fit   | dense-transf |  sparse-fit  | sparse-transf
-----------------|--------------|--------------|--------------|--------------
Bernouilli       |   6.9329s    |   13.7461s   |   7.0625s    |   0.3768s
Gaussian         |   9.4521s    |  148.2662s   |   9.9962s    |   1.5707s

Still TODO:
    - add args to select either sparse or dense problem.

"""
from __future__ import division
from __future__ import print_function

import gc
import sys
import optparse
from datetime import datetime
import collections

import numpy as np
import scipy.sparse as sp

from sklearn import clone
from sklearn.random_projection import (
                                        BernouilliRandomProjection,
                                        GaussianRandomProjection,
                                       )


def type_auto_or_float(val):
    if val == "auto":
        return "auto"
    else:
        return float(val)


def type_auto_or_int(val):
    if val == "auto":
        return "auto"
    else:
        return int(val)


def compute_time(t_start, delta):
    mu_second = 0.0 + 10 ** 6  # number of microseconds in a second

    return delta.seconds + delta.microseconds / mu_second


def bench_scikit_transformer(X, transfomer):
    gc.collect()

    clf = clone(transfomer)

    # start time
    t_start = datetime.now()
    clf.fit(X)
    delta = (datetime.now() - t_start)
    # stop time
    time_to_fit = compute_time(t_start, delta)

    # start time
    t_start = datetime.now()
    clf.transform(X)
    delta = (datetime.now() - t_start)
    # stop time
    time_to_transform = compute_time(t_start, delta)

    return time_to_fit, time_to_transform


# Make some random data with uniformly located non zero entries with
# gaussian distributed values
def make_sparse_random_data(n_samples, n_features, n_nonzeros,
                            random_state=None):
    rng = np.random.RandomState(random_state)
    data_coo = sp.coo_matrix(
        (rng.randn(n_nonzeros),
        (rng.randint(n_samples, size=n_nonzeros),
         rng.randint(n_features, size=n_nonzeros))),
        shape=(n_samples, n_features))
    return data_coo.toarray(), data_coo.tocsr()


def print_row(clf_type, time_dense_fit, time_dense_transform, time_sparse_fit,
              time_sparse_transform):
    print("%s \t | %s | %s | %s | %s" % (clf_type.ljust(12),
                              ("%.4fs" % time_dense_fit).center(12),
                              ("%.4fs" % time_dense_transform).center(12),
                              ("%.4fs" % time_sparse_fit).center(12),
                              ("%.4fs" % time_sparse_transform).center(12)))

if __name__ == "__main__":
    ###########################################################################
    # Option parser
    ###########################################################################
    op = optparse.OptionParser()
    op.add_option("--n-times",
                  dest="n_times", default=10, type=int,
                  help="Bench results are average over n_times experiments")

    op.add_option("--n-features",
                  dest="n_features", default=5 * 10 ** 4, type=int,
                  help="Number of features in the benchmarks")

    op.add_option("--n-components",
                  dest="n_components", default=10 ** 3,
                  help="Size of the random subspace."
                       "('auto' or int > 0)")

    op.add_option("--ratio-nonzeros",
                  dest="ratio_nonzeros", default=10 ** -3, type=float,
                  help="Number of features in the benchmarks")

    op.add_option("--n-samples",
                  dest="n_samples", default=1000, type=int,
                  help="Number of samples in the benchmarks")

    op.add_option("--random-seed",
                  dest="random_seed", default=13, type=int,
                  help="Seed used by the random number generators.")

    op.add_option("--density",
                  dest="density", default=1 / 3,
                  help="Density used by the sparse random projection."
                       "('auto' or float (0.0, 1.0]")

    op.add_option("--eps",
                  dest="eps", default=0.1, type=float,
                  help="See the documentation of the underlying transformers.")

    op.add_option("--transformers",
                  dest="selected_transformers",
                  default='Gaussian,Bernouilli',
                  type=str,
                  help="Comma-separated list of transformer to benchmark. "
                       "Default: %default. Available: "
                       "GaussianRandomProjection,BernouilliRandomProjection")

    (opts, args) = op.parse_args()
    if len(args) > 0:
        op.error("this script takes no arguments.")
        sys.exit(1)
    opts.n_components = type_auto_or_int(opts.n_components)
    opts.density = type_auto_or_float(opts.density)
    selected_transformers = opts.selected_transformers.split(',')

    ###########################################################################
    # Generate dataset
    ###########################################################################
    n_nonzeros = int(opts.ratio_nonzeros * opts.n_features)

    print('=' * 80)
    print('Datasets statics')
    print('=' * 80)
    print('n_features \t= %s' % opts.n_features)
    print('n_samples \t= %s' % opts.n_samples)
    print('n_elements \t= %s' % (opts.n_features * opts.n_samples))
    print('n_nonzeros \t= %s per feature' % n_nonzeros)
    print('ratio_nonzeros \t= %s' % opts.ratio_nonzeros)
    print('n_components \t= %s' % opts.n_components)

    print("")
    print("Generate dataset benchmarks...")
    X_dense, X_sparse = make_sparse_random_data(opts.n_samples,
                                                opts.n_features,
                                                n_nonzeros,
                                                random_state=opts.random_seed)
    print("")

    ###########################################################################
    # Set transfomer input
    ###########################################################################
    transformers = {}

    ###########################################################################
    # Set GaussianRandomProjection input
    gaussian_params = {
        "n_components": opts.n_components,
        "random_state": opts.random_seed
    }
    transformers["Gaussian"] = \
        GaussianRandomProjection(** gaussian_params)

    ###########################################################################
    # Set BernouilliRandomProjection input
    bernouilli_params = {
        "n_components": opts.n_components,
        "random_state": opts.random_seed,
        "density": opts.density,
        "density": opts.eps,
    }

    transformers["Bernouilli"] = \
        BernouilliRandomProjection(** bernouilli_params)

    ###########################################################################
    # Perform benchmark
    ###########################################################################
    time_dense_fit = collections.defaultdict(list)
    time_dense_transform = collections.defaultdict(list)
    time_sparse_fit = collections.defaultdict(list)
    time_sparse_transform = collections.defaultdict(list)

    print('=' * 80)
    print('Benchmarks')
    print('=' * 80)

    for name in selected_transformers:
        print("Perform benchmarks for %s..." % name)

        # Benchmark for dense matrix
        for iteration in xrange(opts.n_times):
            print("iter dense %s" % iteration)
            time_fit, time_transform = bench_scikit_transformer(X_dense,
                transformers[name])
            time_dense_fit[name].append(time_fit)
            time_dense_transform[name].append(time_transform)

        # Benchmark for sparse matrix
        for iteration in xrange(opts.n_times):
            print("iter sparse %s" % iteration)
            time_fit, time_transform = bench_scikit_transformer(X_sparse,
                transformers[name])
            time_sparse_fit[name].append(time_fit)
            time_sparse_transform[name].append(time_transform)

    ###########################################################################
    # Print results
    ###########################################################################
    print("")
    print("Transformer performance:")
    print("===========================")
    print("Results are averaged over %s repetition(s)." % opts.n_times)
    print("")
    print("%s \t | %s | %s | %s | %s" % ("Transformer".ljust(12),
                                         "dense-fit".center(12),
                                         "dense-transf".center(12),
                                         "sparse-fit".center(12),
                                         "sparse-transf".center(12)))
    print(17 * "-" + ("|" + "-" * 14) * 4)

    for name in sorted(selected_transformers):
        print_row(name,
                  np.mean(time_dense_fit[name]),
                  np.mean(time_dense_transform[name]),
                  np.mean(time_sparse_fit[name]),
                  np.mean(time_sparse_transform[name]))

    print("")
    print("")
