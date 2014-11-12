"""
Tests for the birch clustering algorithm.
"""

from scipy import sparse
import numpy as np
import warnings

from .common import generate_clustered_data
from sklearn.cluster.birch import Birch
from sklearn.cluster.hierarchical import AgglomerativeClustering
from sklearn.datasets import make_blobs
from sklearn.linear_model import ElasticNet
from sklearn.metrics import pairwise_distances_argmin

from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns


def test_n_samples_leaves_roots():
    """Sanity check for the number of samples in leaves and roots"""
    X, y = make_blobs(n_samples=10)
    brc = Birch()
    brc.fit(X)
    n_samples_root = sum([sc.n_ for sc in brc.root_.subclusters_])
    n_samples_leaves = sum([sc.n_ for leaf in brc.get_leaves()
                            for sc in leaf.subclusters_])
    assert_equal(n_samples_leaves, X.shape[0])
    assert_equal(n_samples_root, X.shape[0])


def test_partial_fit():
    """Test that fit is equivalent to calling partial_fit multiple times""" 
    X, y = make_blobs(n_samples=100)
    brc = Birch(n_clusters=3)
    brc.fit(X)
    brc_partial = Birch()
    brc_partial.partial_fit(X[:50])
    brc_partial.partial_fit(X[50:])
    assert_array_equal(brc_partial.cluster_centers_, brc.cluster_centers_)
    assert_equal(len(brc.cluster_centers_), 3)


def test_birch_predict():
    """Test the predict method predicts the nearest centroid."""
    rng = np.random.RandomState(0)
    X = generate_clustered_data(n_clusters=3, n_features=3,
                                n_samples_per_cluster=10)

    # n_samples * n_samples_per_cluster
    shuffle_indices = np.arange(30)
    rng.shuffle(shuffle_indices)
    X_shuffle = X[shuffle_indices, :]
    brc = Birch(n_clusters=4, threshold=1.)
    brc.fit(X_shuffle)
    centroids = brc.cluster_centers_
    assert_array_equal(brc.labels_, brc.predict(X_shuffle))
    nearest_centroid = pairwise_distances_argmin(X_shuffle, centroids)
    assert_array_equal(nearest_centroid, brc.labels_)


def test_n_clusters():
    """Test that n_clusters param works properly"""
    X, y = make_blobs(n_samples=100, centers=10)
    brc1 = Birch(n_clusters=10)
    brc1.fit(X)
    assert_equal(len(brc1.cluster_centers_), 10)

    # Test that n_clusters = Agglomerative Clustering gives
    # the same results.
    gc = AgglomerativeClustering(n_clusters=10)
    brc2 = Birch(n_clusters=gc)
    brc2.fit(X)
    assert_equal(len(brc2.cluster_centers_), 10)
    assert_array_equal(brc1.cluster_centers_, brc2.cluster_centers_)
    assert_array_equal(brc1.labels_, brc2.labels_)

    # Test that the wrong global clustering step raises an Error.
    clf = ElasticNet()
    brc3 = Birch(n_clusters=clf)
    assert_raises(ValueError, brc3.fit, X)

    # Test that a small number of clusters raises a warning.
    brc4 = Birch(threshold=10000.)
    assert_warns(UserWarning, brc4.fit, X)


def test_sparse_X():
    """Test that sparse and dense data give same results"""
    X, y = make_blobs(n_samples=100, centers=10)
    brc = Birch(n_clusters=10)
    brc.fit(X)

    csr = sparse.csr_matrix(X)
    brc_sparse = Birch(n_clusters=10)
    brc_sparse.fit(X)

    assert_array_equal(brc.labels_, brc_sparse.labels_)
    assert_array_equal(brc.cluster_centers_, brc_sparse.cluster_centers_)


def check_branching_factor(node, branching_factor):
    subclusters = node.subclusters_
    if node.is_leaf:
        assert_greater_equal(branching_factor, len(subclusters))
        return
    for cluster in subclusters:
        check_branching_factor(cluster.child_, branching_factor)


def test_branching_factor():
    """Test that nodes have at max branching_factor number of subclusters"""
    X, y = make_blobs()
    branching_factor = 9

    # Purposefully set a low threshold to maximize the subclusters.
    brc = Birch(n_clusters=None, branching_factor=branching_factor,
                threshold=0.01)
    brc.fit(X)
    check_branching_factor(brc.root_, branching_factor)

    # Purposefully set a low threshold to maximize the subclusters.
    brc = Birch(n_clusters=3, branching_factor=branching_factor,
                threshold=0.01)
    brc.fit(X)
    check_branching_factor(brc.root_, branching_factor)
