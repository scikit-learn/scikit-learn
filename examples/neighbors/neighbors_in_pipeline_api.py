"""
=====================================
Pipelining nearest neighbors and TSNE
=====================================

This examples demonstrates how to precompute the k nearest neighbors before
using them in TSNE. TSNE can compute the nearest neighbors internally, but
precomputing them can have several benefits, such as finer control, caching for
multiple use, or custom implementations.

This example presents how to chain KNeighborsTransformer and TSNE in a
pipeline, and how to wrap the package `annoy` to replace KNeighborsTransformer
and perform approximate nearest neighbors. This package can be installed with
`pip install annoy`.

Expected output:

Benchmarking on MNIST_2000:
---------------------------
AnnoyTransformer:                    1.885 sec
KNeighborsTransformer:               1.550 sec
TSNE with AnnoyTransformer:          12.408 sec
TSNE with KNeighborsTransformer:     11.856 sec
TSNE with internal NearestNeighbors: 11.884 sec

Benchmarking on MNIST_10000:
----------------------------
AnnoyTransformer:                    12.613 sec
KNeighborsTransformer:               44.348 sec
TSNE with AnnoyTransformer:          75.256 sec
TSNE with KNeighborsTransformer:     110.472 sec
TSNE with internal NearestNeighbors: 110.527 sec
"""
# Author: Tom Dupre la Tour
#
# License: BSD 3 clause
from __future__ import print_function
import time

import annoy
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


class AnnoyTransformer(BaseEstimator, TransformerMixin):
    """Wrapper for using annoy.AnnoyIndex as sklearn's KNeighborsTransformer"""

    def __init__(self, n_neighbors=5, metric='euclidean', n_trees=10,
                 search_k=-1, include_self=False):
        self.n_neighbors = n_neighbors
        self.n_trees = n_trees
        self.search_k = search_k
        self.metric = metric
        self.include_self = include_self

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
            if self.include_self:
                n_neighbors = self.n_neighbors
            else:
                # add a neighbors and remove it later
                n_neighbors = self.n_neighbors + 1

            for i in range(self.annoy_.get_n_items()):
                ind, dist = self.annoy_.get_nns_by_item(
                    i, n_neighbors, self.search_k, include_distances=True)

                if self.include_self:
                    indices[i], distances[i] = ind, dist
                else:
                    indices[i], distances[i] = ind[1:], dist[1:]
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
        self.fit(X)
        return self.transform(X=None)


def test_annoy_transformer():
    """Test that AnnoyTransformer and KNeighborsTransformer give same results
    """
    X = np.random.RandomState(42).randn(10, 2)

    for include_self in [True, False]:
        ann = AnnoyTransformer(include_self=include_self)
        Xt = ann.fit_transform(X)

        knn = KNeighborsTransformer(mode='distance', include_self=include_self)
        Xt2 = knn.fit_transform(X)

        assert_array_almost_equal(Xt.toarray(), Xt2.toarray(), decimal=5)


def load_mnist(n_samples):
    """Load MNIST, shuffle the data, and return only n_samples."""
    mnist = fetch_mldata('MNIST original')
    X, y = shuffle(mnist.data, mnist.target, random_state=42)
    return X[:n_samples], y[:n_samples]


def run_benchmark():
    datasets = {
        'MNIST_2000': load_mnist(n_samples=2000),
        'MNIST_10000': load_mnist(n_samples=10000),
    }

    n_iter = 250
    perplexity = 30
    # TSNE requires a certain number of neighbors which depends on the
    # perplexity parameter
    n_neighbors = int(3. * perplexity + 1)

    transformers = {
        'AnnoyTransformer':  # without TSNE, for speed benchmark only
        AnnoyTransformer(n_neighbors=n_neighbors, include_self=False,
                         metric='euclidean'),

        'KNeighborsTransformer':  # without TSNE, for speed benchmark only
        KNeighborsTransformer(n_neighbors=n_neighbors, mode='distance',
                              metric='sqeuclidean', include_self=False),

        'TSNE with AnnoyTransformer':
        make_pipeline(
            AnnoyTransformer(n_neighbors=n_neighbors, include_self=False,
                             metric='euclidean'),
            TSNE(metric='precomputed', perplexity=perplexity,
                 method="barnes_hut", random_state=42, n_iter=n_iter),
        ),

        'TSNE with KNeighborsTransformer':
        make_pipeline(
            KNeighborsTransformer(n_neighbors=n_neighbors, mode='distance',
                                  metric='sqeuclidean', include_self=False),
            TSNE(metric='precomputed', perplexity=perplexity,
                 method="barnes_hut", random_state=42, n_iter=n_iter),
        ),

        'TSNE with internal NearestNeighbors':
        TSNE(metric='sqeuclidean', perplexity=perplexity, method="barnes_hut",
             random_state=42, n_iter=n_iter),
    }

    # init the plot
    nrows = len(datasets)
    ncols = np.sum([1 for name in transformers.keys() if 'TSNE' in name])
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False,
                             figsize=(5 * ncols, 4 * nrows))
    axes = axes.ravel()
    i_ax = 0

    for dataset_name, (X, y) in datasets.items():

        msg = 'Benchmarking on %s:' % dataset_name
        print('\n%s\n%s' % (msg, '-' * len(msg)))

        for transformer_name, transformer in transformers.items():
            start = time.time()
            Xt = transformer.fit_transform(X)
            duration = time.time() - start

            # print the duration report
            longest = np.max([len(n) for n in transformers.keys()])
            whitespaces = ' ' * (longest - len(transformer_name))
            print('%s: %s%.3f sec' % (transformer_name, whitespaces, duration))

            # plot TSNE embedding which should be very similar across methods
            if 'TSNE' in transformer_name:
                axes[i_ax].set_title(transformer_name + ' on ' + dataset_name)
                axes[i_ax].scatter(Xt[:, 0], Xt[:, 1], c=y,
                                   cmap=plt.cm.viridis)
                axes[i_ax].xaxis.set_major_formatter(NullFormatter())
                axes[i_ax].yaxis.set_major_formatter(NullFormatter())
                axes[i_ax].axis('tight')
                i_ax += 1

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    test_annoy_transformer()

    run_benchmark()
