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


try:
    from memory_profiler import profile
except ImportError:
    def profile(f):
        return f

try:
    # If you want comparison with the reference implementation, build the
    # binary from source (https://github.com/lvdmaaten/bhtsne) in the folder
    # benchmarks/bhtsne and add an empty `__init__.py` file in the floder.
    #
    #  $ git clone git@github.com:lvdmaaten/bhtsne.git
    #  $ cd bhtsne
    #  $ g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2
    #  $ touch __init__.py
    #  $ cd ..
    #
    from bhtsne.bhtsne import run_bh_tsne

    @profile
    def bhtsne(data, **kwargs):
        return run_bh_tsne(data, **kwargs)
except ImportError:
    bhtsne = None


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


@profile
def tsne_fit_transform(model, data):
    return model.fit_transform(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Benchmark for t-SNE')
    parser.add_argument('--order', type=str, default='C',
                        help='Order of the input data')
    parser.add_argument('--perplexity', type=float, default=30)
    parser.add_argument('--log-max-nsamples', type=int, default=5)
    parser.add_argument('--verbose', type=int, default=2)
    args = parser.parse_args()

    X, y = load_data(order=args.order)

    results = []
    basename, _ = os.path.splitext(__file__)
    log_filename = basename + '.json'
    for n in [100, 1000, 5000, 10000]:
        X_train = X[:n]
        n = X_train.shape[0]
        print("Fitting TSNE on %d samples..." % n)
        tsne = TSNE(n_components=2, init='pca', perplexity=args.perplexity,
                    verbose=0)
        t0 = time()
        X_embedded = tsne_fit_transform(tsne, X_train)
        duration = time() - t0
        print("Fitting T-SNE on %d samples took %0.3fs" % (n, duration))
        results.append(dict(method="TSNE", duration=duration, n_samples=n))
        with open(log_filename, 'w', encoding='utf-8') as f:
            json.dump(results, f)
        np.save('mnist_tsne_%d.npy' % n, X_embedded)

        if bhtsne is not None:
            t0 = time()
            X_embedded = bhtsne(X_train, initial_dims=X_train.shape[1],
                                perplexity=args.perplexity, verbose=False)
            duration = time() - t0
            print("Fitting bhtsne on %d samples took %0.3fs" % (n, duration))
            results.append(dict(method="bhtsne", duration=duration,
                                n_samples=n))
            with open(log_filename, 'w', encoding='utf-8') as f:
                json.dump(results, f)
            np.save('mnist_bhtsne_%d.npy' % n, X_embedded)
