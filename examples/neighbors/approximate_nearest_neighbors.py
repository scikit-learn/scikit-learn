"""
=====================================
Approximate nearest neighbors in TSNE
=====================================

This example presents how to chain KNeighborsTransformer and TSNE in a pipeline.
It also shows how to wrap the packages `nmslib` and `pynndescent` to replace
KNeighborsTransformer and perform approximate nearest neighbors. These packages
can be installed with `pip install nmslib pynndescent`.

Note: In KNeighborsTransformer we use the definition which includes each
training point as its own neighbor in the count of `n_neighbors`, and for
compatibility reasons, one extra neighbor is computed when `mode == 'distance'`.
Please note that we do the same in the proposed wrappers.

Sample output::

    Benchmarking on MNIST_2000:
    ---------------------------
    NMSlibTransformer:                   0.144 sec
    KNeighborsTransformer:               0.090 sec
    PyNNDescentTransformer:              23.402 sec
    TSNE with NMSlibTransformer:         2.592 sec
    TSNE with KNeighborsTransformer:     2.338 sec
    TSNE with PyNNDescentTransformer:    6.288 sec
    TSNE with internal NearestNeighbors: 2.364 sec

    Benchmarking on MNIST_10000:
    ----------------------------
    NMSlibTransformer:                   1.098 sec
    KNeighborsTransformer:               1.264 sec
    PyNNDescentTransformer:              7.170 sec
    TSNE with NMSlibTransformer:         15.281 sec
    TSNE with KNeighborsTransformer:     15.400 sec
    TSNE with PyNNDescentTransformer:    28.782 sec
    TSNE with internal NearestNeighbors: 15.573 sec


Note that the prediction speed KNeighborsTransformer was optimized in
scikit-learn 1.1 and therefore approximate methods are not necessarily faster
because computing the index takes time and can nullify the gains obtained at
prediction time.

"""

# Author: Tom Dupre la Tour
#
# License: BSD 3 clause
import joblib
import time
import sys

try:
    import nmslib
except ImportError:
    print("The package 'nmslib' is required to run this example.")
    sys.exit()

try:
    from pynndescent import PyNNDescentTransformer
except ImportError:
    print("The package 'pynndescent' is required to run this example.")
    sys.exit()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.sparse import csr_matrix

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.neighbors import KNeighborsTransformer
from sklearn.datasets import fetch_openml
from sklearn.pipeline import make_pipeline
from sklearn.manifold import TSNE
from sklearn.utils import shuffle


class NMSlibTransformer(TransformerMixin, BaseEstimator):
    """Wrapper for using nmslib as sklearn's KNeighborsTransformer"""

    def __init__(self, n_neighbors=5, metric="euclidean", method="sw-graph", n_jobs=-1):
        self.n_neighbors = n_neighbors
        self.method = method
        self.metric = metric
        self.n_jobs = n_jobs

    def fit(self, X):
        self.n_samples_fit_ = X.shape[0]

        # see more metric in the manual
        # https://github.com/nmslib/nmslib/tree/master/manual
        space = {
            "euclidean": "l2",
            "cosine": "cosinesimil",
            "l1": "l1",
            "l2": "l2",
        }[self.metric]

        self.nmslib_ = nmslib.init(method=self.method, space=space)
        self.nmslib_.addDataPointBatch(X.copy())
        self.nmslib_.createIndex()
        return self

    def transform(self, X):
        n_samples_transform = X.shape[0]

        # For compatibility reasons, as each sample is considered as its own
        # neighbor, one extra neighbor will be computed.
        n_neighbors = self.n_neighbors + 1

        if self.n_jobs < 0:
            # Same handling as done in joblib for negative values of n_jobs:
            # in particular, `n_jobs == -1` means "as many threads as CPUs".
            num_threads = joblib.cpu_count() + self.n_jobs + 1
        else:
            num_threads = self.n_jobs

        results = self.nmslib_.knnQueryBatch(
            X.copy(), k=n_neighbors, num_threads=num_threads
        )
        indices, distances = zip(*results)
        indices, distances = np.vstack(indices), np.vstack(distances)

        indptr = np.arange(0, n_samples_transform * n_neighbors + 1, n_neighbors)
        kneighbors_graph = csr_matrix(
            (distances.ravel(), indices.ravel(), indptr),
            shape=(n_samples_transform, self.n_samples_fit_),
        )

        return kneighbors_graph


def load_mnist(n_samples):
    """Load MNIST, shuffle the data, and return only n_samples."""
    mnist = fetch_openml("mnist_784", as_frame=False, parser="pandas")
    X, y = shuffle(mnist.data, mnist.target, random_state=2)
    return X[:n_samples] / 255, y[:n_samples]


def run_benchmark():
    datasets = [
        ("MNIST_2000", load_mnist(n_samples=2000)),
        ("MNIST_10000", load_mnist(n_samples=10000)),
    ]

    n_iter = 500
    perplexity = 30
    metric = "euclidean"
    # TSNE requires a certain number of neighbors which depends on the
    # perplexity parameter.
    # Add one since we include each sample as its own neighbor.
    n_neighbors = int(3.0 * perplexity + 1) + 1

    tsne_params = dict(
        init="random",  # pca not supported for sparse matrices
        perplexity=perplexity,
        method="barnes_hut",
        random_state=42,
        n_iter=n_iter,
        learning_rate="auto",
    )

    transformers = [
        (
            "NMSlibTransformer",
            NMSlibTransformer(n_neighbors=n_neighbors, metric=metric),
        ),
        (
            "KNeighborsTransformer",
            KNeighborsTransformer(
                n_neighbors=n_neighbors, mode="distance", metric=metric
            ),
        ),
        (
            "PyNNDescentTransformer",
            PyNNDescentTransformer(n_neighbors=n_neighbors, metric=metric),
        ),
        (
            "TSNE with NMSlibTransformer",
            make_pipeline(
                NMSlibTransformer(n_neighbors=n_neighbors, metric=metric),
                TSNE(metric="precomputed", **tsne_params),
            ),
        ),
        (
            "TSNE with KNeighborsTransformer",
            make_pipeline(
                KNeighborsTransformer(
                    n_neighbors=n_neighbors, mode="distance", metric=metric
                ),
                TSNE(metric="precomputed", **tsne_params),
            ),
        ),
        (
            "TSNE with PyNNDescentTransformer",
            make_pipeline(
                PyNNDescentTransformer(n_neighbors=n_neighbors, metric=metric),
                TSNE(metric="precomputed", **tsne_params),
            ),
        ),
        ("TSNE with internal NearestNeighbors", TSNE(metric=metric, **tsne_params)),
    ]

    # init the plot
    nrows = len(datasets)
    ncols = np.sum([1 for name, model in transformers if "TSNE" in name])
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, squeeze=False, figsize=(5 * ncols, 4 * nrows)
    )
    axes = axes.ravel()
    i_ax = 0

    for dataset_name, (X, y) in datasets:

        msg = "Benchmarking on %s:" % dataset_name
        print("\n%s\n%s" % (msg, "-" * len(msg)))

        for transformer_name, transformer in transformers:
            start = time.time()
            Xt = transformer.fit_transform(X)
            duration = time.time() - start

            # print the duration report
            longest = np.max([len(name) for name, model in transformers])
            whitespaces = " " * (longest - len(transformer_name))
            print("%s: %s%.3f sec" % (transformer_name, whitespaces, duration))

            # plot TSNE embedding which should be very similar across methods
            if "TSNE" in transformer_name:
                axes[i_ax].set_title(transformer_name + "\non " + dataset_name)
                axes[i_ax].scatter(
                    Xt[:, 0],
                    Xt[:, 1],
                    c=y.astype(np.int32),
                    alpha=0.2,
                    cmap=plt.cm.viridis,
                )
                axes[i_ax].xaxis.set_major_formatter(NullFormatter())
                axes[i_ax].yaxis.set_major_formatter(NullFormatter())
                axes[i_ax].axis("tight")
                i_ax += 1

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_benchmark()
