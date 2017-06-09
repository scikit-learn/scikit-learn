"""
=============================
MNIST dataset T-SNE benchmark
=============================

"""
from __future__ import division, print_function

# License: BSD 3 clause

import os
from time import time
import numpy as np
import json
import argparse

from sklearn.externals.joblib import Memory
from sklearn.datasets import fetch_mldata
from sklearn.manifold import TSNE
from sklearn.utils import check_array


memory = Memory('mnist_tsne_benchmark_data', mmap_mode='r')


@memory.cache
def load_data(dtype=np.float32, order='C'):
    """Load the data, then cache and memmap the train/test split"""
    print("Loading dataset...")
    data = fetch_mldata('MNIST original')
    X = check_array(data['data'], dtype=dtype, order=order)
    y = data["target"]

    # Normalize features
    X /= 255
    return X, y


def tsne_fit_transform(model, data):
    return model.fit_transform(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Benchmark for t-SNE')
    parser.add_argument('--order', type=str, default='C',
                        help='Order of the input data')
    parser.add_argument('--perplexity', type=float, default=30)
    parser.add_argument('--bhtsne', action='store_true',
                        help="if set and the reference bhtsne code is "
                        "correctly installed, run it in the benchmark.")
    parser.add_argument('--all', action='store_true',
                        help="if set, run the benchmark with the whole MNIST."
                        "dataset. Note that it will take up to 1hour.")
    parser.add_argument('--profile', action='store_true',
                        help="if set, run the benchmark with a memory "
                        "profiler.")
    parser.add_argument('--verbose', type=int, default=0)
    parser.add_argument('--n_jobs', type=int, nargs="+", default=1,
                        help="Number of CPU used to fit sklearn.TSNE")
    args = parser.parse_args()

    X, y = load_data(order=args.order)

    methods = []

    # Put TSNE in methods
    if isinstance(args.n_jobs, int):
        tsne = TSNE(n_components=2, init='pca', perplexity=args.perplexity,
                    verbose=args.verbose, n_jobs=args.n_jobs)
        methods += [("sklearn.TSNE", tsne.fit_transform)]
    elif isinstance(args.n_jobs, list):
        for n_jobs in args.n_jobs:
            tsne = TSNE(n_components=2, init='pca', perplexity=args.perplexity,
                        verbose=args.verbose, n_jobs=n_jobs)
            methods += [("sklearn.TSNE_n_jobs{}".format(n_jobs),
                        tsne.fit_transform)]

    if args.bhtsne:
        try:
            from bhtsne.bhtsne import run_bh_tsne
        except ImportError:
            raise ImportError(
                """
        If you want comparison with the reference implementation, build the
        binary from source (https://github.com/lvdmaaten/bhtsne) in the folder
        benchmarks/bhtsne and add an empty `__init__.py` file in the folder:

        $ git clone git@github.com:lvdmaaten/bhtsne.git
        $ cd bhtsne
        $ g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2
        $ touch __init__.py
        $ cd ..
                   """
            )

        def bhtsne(X):
            """wrapper for LvdM bhtsne implementation."""
            return run_bh_tsne(X, initial_dims=X.shape[1],
                               perplexity=args.perplexity, verbose=False)
        methods += [("LvdM.bhtsne", bhtsne)]

    if args.profile:

        try:
            from memory_profiler import profile
        except ImportError:
            raise ImportError("To run the benchmark with `--profile`, you "
                              "need to install `memory_profiler`. Please "
                              "run `pip install memory_profiler`.")
        methods = [(n, profile(m)) for n, m in methods]

    data_size = [100, 1000, 5000, 10000]
    if args.all:
        data_size += [70000]

    results = []
    basename, _ = os.path.splitext(__file__)
    log_filename = basename + '.json'
    for n in data_size:
        X_train = X[:n]
        n = X_train.shape[0]
        for name, method in methods:
            print("Fitting {} on {} samples...".format(name, n))
            t0 = time()
            X_embedded = method(X_train)
            duration = time() - t0
            print("Fitting {} on {} samples took {:.3f}s"
                  .format(name, n, duration))
            results.append(dict(method=name, duration=duration, n_samples=n))
            with open(log_filename, 'w', encoding='utf-8') as f:
                json.dump(results, f)
            np.save('mnist_{}_{}.npy'.format(name, n), X_embedded)
