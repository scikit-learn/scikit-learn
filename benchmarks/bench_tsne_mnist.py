"""
=============================
MNIST dataset T-SNE benchmark
=============================

"""

# License: BSD 3 clause

import os
import os.path as op
from time import time
import numpy as np
import json
import argparse
from joblib import Memory

from sklearn.datasets import fetch_openml
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.utils import check_array
from sklearn.utils import shuffle as _shuffle
from sklearn.utils._openmp_helpers import _openmp_effective_n_threads

LOG_DIR = "mnist_tsne_output"
if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)


memory = Memory(os.path.join(LOG_DIR, "mnist_tsne_benchmark_data"), mmap_mode="r")


@memory.cache
def load_data(dtype=np.float32, order="C", shuffle=True, seed=0):
    """Load the data, then cache and memmap the train/test split"""
    print("Loading dataset...")
    data = fetch_openml("mnist_784", as_frame=True, parser="pandas")

    X = check_array(data["data"], dtype=dtype, order=order)
    y = data["target"]

    if shuffle:
        X, y = _shuffle(X, y, random_state=seed)

    # Normalize features
    X /= 255
    return X, y


def nn_accuracy(X, X_embedded, k=1):
    """Accuracy of the first nearest neighbor"""
    knn = NearestNeighbors(n_neighbors=1, n_jobs=-1)
    _, neighbors_X = knn.fit(X).kneighbors()
    _, neighbors_X_embedded = knn.fit(X_embedded).kneighbors()
    return np.mean(neighbors_X == neighbors_X_embedded)


def tsne_fit_transform(model, data):
    transformed = model.fit_transform(data)
    return transformed, model.n_iter_


def sanitize(filename):
    return filename.replace("/", "-").replace(" ", "_")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Benchmark for t-SNE")
    parser.add_argument(
        "--order", type=str, default="C", help="Order of the input data"
    )
    parser.add_argument("--perplexity", type=float, default=30)
    parser.add_argument(
        "--bhtsne",
        action="store_true",
        help=(
            "if set and the reference bhtsne code is "
            "correctly installed, run it in the benchmark."
        ),
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help=(
            "if set, run the benchmark with the whole MNIST."
            "dataset. Note that it will take up to 1 hour."
        ),
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="if set, run the benchmark with a memory profiler.",
    )
    parser.add_argument("--verbose", type=int, default=0)
    parser.add_argument(
        "--pca-components",
        type=int,
        default=50,
        help="Number of principal components for preprocessing.",
    )
    args = parser.parse_args()

    print("Used number of threads: {}".format(_openmp_effective_n_threads()))
    X, y = load_data(order=args.order)

    if args.pca_components > 0:
        t0 = time()
        X = PCA(n_components=args.pca_components).fit_transform(X)
        print(
            "PCA preprocessing down to {} dimensions took {:0.3f}s".format(
                args.pca_components, time() - t0
            )
        )

    methods = []

    # Put TSNE in methods
    tsne = TSNE(
        n_components=2,
        init="pca",
        perplexity=args.perplexity,
        verbose=args.verbose,
        n_iter=1000,
    )
    methods.append(("sklearn TSNE", lambda data: tsne_fit_transform(tsne, data)))

    if args.bhtsne:
        try:
            from bhtsne.bhtsne import run_bh_tsne
        except ImportError as e:
            raise ImportError(
                """\
If you want comparison with the reference implementation, build the
binary from source (https://github.com/lvdmaaten/bhtsne) in the folder
benchmarks/bhtsne and add an empty `__init__.py` file in the folder:

$ git clone git@github.com:lvdmaaten/bhtsne.git
$ cd bhtsne
$ g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2
$ touch __init__.py
$ cd ..
"""
            ) from e

        def bhtsne(X):
            """Wrapper for the reference lvdmaaten/bhtsne implementation."""
            # PCA preprocessing is done elsewhere in the benchmark script
            n_iter = -1  # TODO find a way to report the number of iterations
            return (
                run_bh_tsne(
                    X,
                    use_pca=False,
                    perplexity=args.perplexity,
                    verbose=args.verbose > 0,
                ),
                n_iter,
            )

        methods.append(("lvdmaaten/bhtsne", bhtsne))

    if args.profile:

        try:
            from memory_profiler import profile
        except ImportError as e:
            raise ImportError(
                "To run the benchmark with `--profile`, you "
                "need to install `memory_profiler`. Please "
                "run `pip install memory_profiler`."
            ) from e
        methods = [(n, profile(m)) for n, m in methods]

    data_size = [100, 500, 1000, 5000, 10000]
    if args.all:
        data_size.append(70000)

    results = []
    basename = os.path.basename(os.path.splitext(__file__)[0])
    log_filename = os.path.join(LOG_DIR, basename + ".json")
    for n in data_size:
        X_train = X[:n]
        y_train = y[:n]
        n = X_train.shape[0]
        for name, method in methods:
            print("Fitting {} on {} samples...".format(name, n))
            t0 = time()
            np.save(
                os.path.join(LOG_DIR, "mnist_{}_{}.npy".format("original", n)), X_train
            )
            np.save(
                os.path.join(LOG_DIR, "mnist_{}_{}.npy".format("original_labels", n)),
                y_train,
            )
            X_embedded, n_iter = method(X_train)
            duration = time() - t0
            precision_5 = nn_accuracy(X_train, X_embedded)
            print(
                "Fitting {} on {} samples took {:.3f}s in {:d} iterations, "
                "nn accuracy: {:0.3f}".format(name, n, duration, n_iter, precision_5)
            )
            results.append(dict(method=name, duration=duration, n_samples=n))
            with open(log_filename, "w", encoding="utf-8") as f:
                json.dump(results, f)
            method_name = sanitize(name)
            np.save(
                op.join(LOG_DIR, "mnist_{}_{}.npy".format(method_name, n)), X_embedded
            )
