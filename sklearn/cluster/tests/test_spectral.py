"""Testing for Spectral Clustering methods"""

from cPickle import dumps, loads
import nose

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_raises
from scipy import sparse
from sklearn.datasets.samples_generator import make_blobs
from sklearn.utils.testing import assert_greater

from sklearn.cluster import SpectralClustering, spectral_clustering
from sklearn.cluster.spectral import spectral_embedding
from sklearn.metrics import pairwise_distances, adjusted_rand_score


def test_spectral_clustering():
    S = np.array([[1, 5, 2, 1, 0, 0, 0],
                  [5, 1, 3, 1, 0, 0, 0],
                  [2, 3, 1, 1, 0, 0, 0],
                  [1, 1, 1, 1, 2, 1, 1],
                  [0, 0, 0, 2, 2, 3, 2],
                  [0, 0, 0, 1, 3, 1, 4],
                  [0, 0, 0, 1, 2, 4, 1],
                 ])

    for mode in ('arpack', 'lobpcg'):
        for mat in (S, sparse.csr_matrix(S)):
            model = SpectralClustering(random_state=0, n_clusters=2,
                    affinity='precomputed', mode=mode).fit(mat)
            labels = model.labels_
            if labels[0] == 0:
                labels = 1 - labels

            assert_equal(labels, [1, 1, 1, 0, 0, 0, 0])

            model_copy = loads(dumps(model))
            assert_equal(model_copy.n_clusters, model.n_clusters)
            assert_equal(model_copy.mode, model.mode)
            assert_equal(model_copy.random_state.get_state(),
                        model.random_state.get_state())
            assert_equal(model_copy.labels_, model.labels_)


def test_spectral_amg_mode():
    # Test the amg mode of SpectralClustering
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X)  # Distance matrix
    S = np.max(D) - D  # Similarity matrix
    S = sparse.coo_matrix(S)
    try:
        from pyamg import smoothed_aggregation_solver
        amg_loaded = True
    except ImportError:
        amg_loaded = False
    if amg_loaded:
        labels = spectral_clustering(S, n_clusters=len(centers),
                                     random_state=0, mode="amg")
        # We don't care too much that it's good, just that it *worked*.
        # There does have to be some lower limit on the performance though.
        assert_greater(np.mean(labels == true_labels), .3)
    else:
        assert_raises(ValueError, spectral_embedding, S,
                      n_components=len(centers), random_state=0, mode="amg")


def test_spectral_unknown_mode():
    # Test that SpectralClustering fails with an unknown mode set.
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X)  # Distance matrix
    S = np.max(D) - D  # Similarity matrix
    S = sparse.coo_matrix(S)
    assert_raises(ValueError, spectral_clustering, S, n_clusters=2,
                  random_state=0, mode="<unknown>")


def test_spectral_clustering_sparse():
    # We need a large matrice, or the lobpcg solver will fallback to its
    # non-sparse and buggy mode
    S = np.array([[1, 5, 2, 2, 1, 0, 0, 0, 0, 0],
                  [5, 1, 3, 2, 1, 0, 0, 0, 0, 0],
                  [2, 3, 1, 1, 1, 0, 0, 0, 0, 0],
                  [2, 2, 1, 1, 1, 0, 0, 0, 0, 0],
                  [1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
                  [0, 0, 0, 0, 1, 2, 2, 3, 3, 2],
                  [0, 0, 0, 0, 2, 2, 3, 3, 3, 4],
                  [0, 0, 0, 0, 1, 3, 3, 1, 2, 4],
                  [0, 0, 0, 0, 1, 3, 3, 2, 1, 4],
                  [0, 0, 0, 0, 1, 2, 4, 4, 4, 1],
                 ])

    S = sparse.coo_matrix(S)

    labels = SpectralClustering(random_state=0, n_clusters=2,
            affinity='precomputed').fit(S).labels_
    if labels[0] == 0:
        labels = 1 - labels

    assert_greater(np.mean(labels == [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]), .9)


def test_affinities():
    X, y = make_blobs(n_samples=40, random_state=1, centers=[[1, 1], [-1, -1]],
            cluster_std=0.4)
    # nearest neighbors affinity
    sp = SpectralClustering(n_clusters=2, affinity='nearest_neighbors',
            random_state=0)
    labels = sp.fit(X).labels_
    assert_equal(adjusted_rand_score(y, labels), 1)

    sp = SpectralClustering(n_clusters=2, gamma=2, random_state=0)
    labels = sp.fit(X).labels_
    assert_equal(adjusted_rand_score(y, labels), 1)
