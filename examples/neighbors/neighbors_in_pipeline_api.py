"""
=====================================
Pipelining nearest neighbors and TSNE
=====================================

This examples demonstrates how to precompute the k nearest neighbors before
using them in TSNE. TSNE can compute the nearest neighbors internally, but
precomputing them can have several benefits, such as finer control, caching for
multiple use, or custom implementations.

This example presents how to chain KNeighborsTransformer and TSNE in a
pipeline, and how to wrap the packages `annoy` and `nmslib` to replace
KNeighborsTransformer and perform approximate nearest neighbors.
These package can be installed with `pip install annoy nmslib`.

Expected output:

Benchmarking on MNIST_2000:
---------------------------
AnnoyTransformer:                    0.583 sec
NMSlibTransformer:                   0.321 sec
KNeighborsTransformer:               1.225 sec
TSNE with AnnoyTransformer:          4.903 sec
TSNE with NMSlibTransformer:         5.009 sec
TSNE with KNeighborsTransformer:     6.210 sec
TSNE with internal NearestNeighbors: 6.365 sec

Benchmarking on MNIST_10000:
----------------------------
AnnoyTransformer:                    4.457 sec
NMSlibTransformer:                   2.080 sec
KNeighborsTransformer:               30.680 sec
TSNE with AnnoyTransformer:          30.225 sec
TSNE with NMSlibTransformer:         43.295 sec
TSNE with KNeighborsTransformer:     64.845 sec
TSNE with internal NearestNeighbors: 64.984 sec
"""
# Author: Tom Dupre la Tour
#
# License: BSD 3 clause
from __future__ import print_function
import time
import sys

try:
    import annoy
except ImportError:
    print("The package 'annoy' is required to run this example.")
    sys.exit()

try:
    import nmslib
except ImportError:
    print("The package 'nmslib' is required to run this example.")
    sys.exit()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.sparse import csr_matrix

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.neighbors import KNeighborsTransformer
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.datasets import fetch_mldata
from sklearn.pipeline import make_pipeline
from sklearn.manifold import TSNE
from sklearn.utils import shuffle

print(__doc__)


class NMSlibTransformer(BaseEstimator, TransformerMixin):
    """Wrapper for using nmslib as sklearn's KNeighborsTransformer"""

    def __init__(self, n_neighbors=5, metric='euclidean', method='sw-graph',
                 n_jobs=1):
        self.n_neighbors = n_neighbors
        self.method = method
        self.metric = metric
        self.n_jobs = n_jobs

    def fit(self, X):
        self.n_samples_fit = X.shape[0]

        # see more metric in the manual
        # https://github.com/nmslib/nmslib/blob/master/manual/manual.pdf
        space = {
            'euclidean': 'l2',
            'cosine': 'cosinesimil',
            'l1': 'l1',
            'l2': 'l2'
        }[self.metric]

        self.nmslib_ = nmslib.init(method=self.method, space=space)
        self.nmslib_.addDataPointBatch(X)
        self.nmslib_.createIndex()
        return self

    def transform(self, X):
        n_samples_transform = X.shape[0]

        results = self.nmslib_.knnQueryBatch(X, k=self.n_neighbors,
                                             num_threads=self.n_jobs)
        indices, distances = zip(*results)
        indices, distances = np.vstack(indices), np.vstack(distances)

        indptr = np.arange(0, n_samples_transform * self.n_neighbors + 1,
                           self.n_neighbors)
        kneighbors_graph = csr_matrix((distances.ravel(), indices.ravel(),
                                       indptr), shape=(n_samples_transform,
                                                       self.n_samples_fit))

        return kneighbors_graph

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)


class AnnoyTransformer(BaseEstimator, TransformerMixin):
    """Wrapper for using annoy.AnnoyIndex as sklearn's KNeighborsTransformer"""

    def __init__(self, n_neighbors=5, metric='euclidean', n_trees=10,
                 search_k=-1):
        self.n_neighbors = n_neighbors
        self.n_trees = n_trees
        self.search_k = search_k
        self.metric = metric

    def fit(self, X):
        self.n_samples_fit = X.shape[0]
        self.annoy_ = annoy.AnnoyIndex(X.shape[1], metric=self.metric)
        for i, x in enumerate(X):
            self.annoy_.add_item(i, x.tolist())
        self.annoy_.build(self.n_trees)
        return self

    def transform(self, X):
        if X is None:
            n_samples_transform = self.n_samples_fit
        else:
            n_samples_transform = X.shape[0]

        indices = np.empty((n_samples_transform, self.n_neighbors),
                           dtype=np.int)
        distances = np.empty((n_samples_transform, self.n_neighbors))

        if X is None:
            n_neighbors = self.n_neighbors

            for i in range(self.annoy_.get_n_items()):
                ind, dist = self.annoy_.get_nns_by_item(
                    i, n_neighbors, self.search_k, include_distances=True)

                indices[i], distances[i] = ind, dist
        else:
            for i, x in enumerate(X):
                indices[i], distances[i] = self.annoy_.get_nns_by_vector(
                    x.tolist(), self.n_neighbors, self.search_k,
                    include_distances=True)

        indptr = np.arange(0, n_samples_transform * self.n_neighbors + 1,
                           self.n_neighbors)
        kneighbors_graph = csr_matrix((distances.ravel(), indices.ravel(),
                                       indptr), shape=(n_samples_transform,
                                                       self.n_samples_fit))

        return kneighbors_graph

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X=None)


def test_transformers():
    """Test that AnnoyTransformer and KNeighborsTransformer give same results
    """
    X = np.random.RandomState(42).randn(10, 2)

    knn = KNeighborsTransformer(mode='distance')
    Xt = knn.fit_transform(X)

    ann = AnnoyTransformer()
    Xt2 = ann.fit_transform(X)

    nms = NMSlibTransformer()
    Xt3 = nms.fit_transform(X)

    assert_array_almost_equal(Xt.toarray(), Xt2.toarray(), decimal=5)
    assert_array_almost_equal(Xt.toarray(), Xt3.toarray(), decimal=5)


def load_mnist(n_samples):
    """Load MNIST, shuffle the data, and return only n_samples."""
    mnist = fetch_mldata('MNIST original')
    X, y = shuffle(mnist.data, mnist.target, random_state=42)
    return X[:n_samples], y[:n_samples]


def run_benchmark():
    datasets = [
        ('MNIST_2000', load_mnist(n_samples=2000)),
        ('MNIST_10000', load_mnist(n_samples=10000)),
    ]

    n_iter = 250
    perplexity = 30
    # TSNE requires a certain number of neighbors which depends on the
    # perplexity parameter
    # Add one since we include each sample as its own neighbor.
    n_neighbors = int(3. * perplexity + 1) + 1

    transformers = [
        ('AnnoyTransformer', AnnoyTransformer(n_neighbors=n_neighbors,
                                              metric='euclidean')),
        ('NMSlibTransformer', NMSlibTransformer(n_neighbors=n_neighbors,
                                                metric='euclidean')),
        ('KNeighborsTransformer', KNeighborsTransformer(
            n_neighbors=n_neighbors, mode='distance', metric='sqeuclidean')),

        ('TSNE with AnnoyTransformer', make_pipeline(
            AnnoyTransformer(n_neighbors=n_neighbors, metric='euclidean'),
            TSNE(metric='precomputed', perplexity=perplexity,
                 method="barnes_hut", random_state=42, n_iter=n_iter), )),
        ('TSNE with NMSlibTransformer', make_pipeline(
            NMSlibTransformer(n_neighbors=n_neighbors, metric='euclidean'),
            TSNE(metric='precomputed', perplexity=perplexity,
                 method="barnes_hut", random_state=42, n_iter=n_iter), )),
        ('TSNE with KNeighborsTransformer', make_pipeline(
            KNeighborsTransformer(n_neighbors=n_neighbors, mode='distance',
                                  metric='sqeuclidean'),
            TSNE(metric='precomputed', perplexity=perplexity,
                 method="barnes_hut", random_state=42, n_iter=n_iter), )),
        ('TSNE with internal NearestNeighbors',
         TSNE(metric='sqeuclidean', perplexity=perplexity, method="barnes_hut",
              random_state=42, n_iter=n_iter)),
    ]

    # init the plot
    nrows = len(datasets)
    ncols = np.sum([1 for name, model in transformers if 'TSNE' in name])
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False,
                             figsize=(5 * ncols, 4 * nrows))
    axes = axes.ravel()
    i_ax = 0

    for dataset_name, (X, y) in datasets:

        msg = 'Benchmarking on %s:' % dataset_name
        print('\n%s\n%s' % (msg, '-' * len(msg)))

        for transformer_name, transformer in transformers:
            start = time.time()
            Xt = transformer.fit_transform(X)
            duration = time.time() - start

            # print the duration report
            longest = np.max([len(name) for name, model in transformers])
            whitespaces = ' ' * (longest - len(transformer_name))
            print('%s: %s%.3f sec' % (transformer_name, whitespaces, duration))

            # plot TSNE embedding which should be very similar across methods
            if 'TSNE' in transformer_name:
                axes[i_ax].set_title(transformer_name + ' on ' + dataset_name)
                axes[i_ax].scatter(Xt[:, 0], Xt[:, 1], c=y, alpha=0.2,
                                   cmap=plt.cm.viridis)
                axes[i_ax].xaxis.set_major_formatter(NullFormatter())
                axes[i_ax].yaxis.set_major_formatter(NullFormatter())
                axes[i_ax].axis('tight')
                i_ax += 1

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    test_transformers()
    run_benchmark()
