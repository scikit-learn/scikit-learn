"""
===========================
Random projection benchmark
===========================

Benchmarks for random projections.

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
from sklearn.random_projection import (SparseRandomProjection,
                                       GaussianRandomProjection,
                                       johnson_lindenstrauss_min_dim)


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
# Gaussian distributed values
def make_sparse_random_data(n_samples, n_features, n_nonzeros,
                            random_state=None):
    rng = np.random.RandomState(random_state)
    data_coo = sp.coo_matrix(
        (rng.randn(n_nonzeros),
        (rng.randint(n_samples, size=n_nonzeros),
         rng.randint(n_features, size=n_nonzeros))),
        shape=(n_samples, n_features))
    return data_coo.toarray(), data_coo.tocsr()


def print_row(clf_type, time_fit, time_transform):
    print("%s | %s | %s" % (clf_type.ljust(30),
                           ("%.4fs" % time_fit).center(12),
                           ("%.4fs" % time_transform).center(12)))


if __name__ == "__main__":
    ###########################################################################
    # Option parser
    ###########################################################################
    op = optparse.OptionParser()
    op.add_option("--n-times",
                  dest="n_times", default=5, type=int,
                  help="Benchmark results are average over n_times experiments")

    op.add_option("--n-features",
                  dest="n_features", default=10 ** 4, type=int,
                  help="Number of features in the benchmarks")

    op.add_option("--n-components",
                  dest="n_components", default="auto",
                  help="Size of the random subspace."
                       " ('auto' or int > 0)")

    op.add_option("--ratio-nonzeros",
                  dest="ratio_nonzeros", default=10 ** -3, type=float,
                  help="Number of features in the benchmarks")

    op.add_option("--n-samples",
                  dest="n_samples", default=500, type=int,
                  help="Number of samples in the benchmarks")

    op.add_option("--random-seed",
                  dest="random_seed", default=13, type=int,
                  help="Seed used by the random number generators.")

    op.add_option("--density",
                  dest="density", default=1 / 3,
                  help="Density used by the sparse random projection."
                       " ('auto' or float (0.0, 1.0]")

    op.add_option("--eps",
                  dest="eps", default=0.5, type=float,
                  help="See the documentation of the underlying transformers.")

    op.add_option("--transformers",
                  dest="selected_transformers",
                  default='GaussianRandomProjection,SparseRandomProjection',
                  type=str,
                  help="Comma-separated list of transformer to benchmark. "
                       "Default: %default. Available: "
                       "GaussianRandomProjection,SparseRandomProjection")

    op.add_option("--dense",
                  dest="dense",
                  default=False,
                  action="store_true",
                  help="Set input space as a dense matrix.")

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

    print('Dataset statics')
    print("===========================")
    print('n_samples \t= %s' % opts.n_samples)
    print('n_features \t= %s' % opts.n_features)
    if opts.n_components == "auto":
        print('n_components \t= %s (auto)' %
              johnson_lindenstrauss_min_dim(n_samples=opts.n_samples,
                                            eps=opts.eps))
    else:
        print('n_components \t= %s' % opts.n_components)
    print('n_elements \t= %s' % (opts.n_features * opts.n_samples))
    print('n_nonzeros \t= %s per feature' % n_nonzeros)
    print('ratio_nonzeros \t= %s' % opts.ratio_nonzeros)
    print('')

    ###########################################################################
    # Set transformer input
    ###########################################################################
    transformers = {}

    ###########################################################################
    # Set GaussianRandomProjection input
    gaussian_matrix_params = {
        "n_components": opts.n_components,
        "random_state": opts.random_seed
    }
    transformers["GaussianRandomProjection"] = \
        GaussianRandomProjection(**gaussian_matrix_params)

    ###########################################################################
    # Set SparseRandomProjection input
    sparse_matrix_params = {
        "n_components": opts.n_components,
        "random_state": opts.random_seed,
        "density": opts.density,
        "eps": opts.eps,
    }

    transformers["SparseRandomProjection"] = \
        SparseRandomProjection(**sparse_matrix_params)

    ###########################################################################
    # Perform benchmark
    ###########################################################################
    time_fit = collections.defaultdict(list)
    time_transform = collections.defaultdict(list)

    print('Benchmarks')
    print("===========================")
    print("Generate dataset benchmarks... ", end="")
    X_dense, X_sparse = make_sparse_random_data(opts.n_samples,
                                                opts.n_features,
                                                n_nonzeros,
                                                random_state=opts.random_seed)
    X = X_dense if opts.dense else X_sparse
    print("done")

    for name in selected_transformers:
        print("Perform benchmarks for %s..." % name)

        for iteration in range(opts.n_times):
            print("\titer %s..." % iteration, end="")
            time_to_fit, time_to_transform = bench_scikit_transformer(X_dense,
              transformers[name])
            time_fit[name].append(time_to_fit)
            time_transform[name].append(time_to_transform)
            print("done")

    print("")

    ###########################################################################
    # Print results
    ###########################################################################
    print("Script arguments")
    print("===========================")
    arguments = vars(opts)
    print("%s \t | %s " % ("Arguments".ljust(16),
                           "Value".center(12),))
    print(25 * "-" + ("|" + "-" * 14) * 1)
    for key, value in arguments.items():
        print("%s \t | %s " % (str(key).ljust(16),
                               str(value).strip().center(12)))
    print("")

    print("Transformer performance:")
    print("===========================")
    print("Results are averaged over %s repetition(s)." % opts.n_times)
    print("")
    print("%s | %s | %s" % ("Transformer".ljust(30),
                            "fit".center(12),
                            "transform".center(12)))
    print(31 * "-" + ("|" + "-" * 14) * 2)

    for name in sorted(selected_transformers):
        print_row(name,
                  np.mean(time_fit[name]),
                  np.mean(time_transform[name]))

    print("")
    print("")
