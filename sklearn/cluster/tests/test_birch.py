"""
Tests for the birch clustering algorithm.
"""

import numpy as np

from sklearn.cluster.birch import Birch
from sklearn.datasets import make_blobs

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal


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
    X, y = make_blobs(n_samples=100)
    brc = Birch()
    brc.fit(X)
    brc_partial = Birch()
    brc_partial.partial_fit(X[:50])
    brc_partial.partial_fit(X[50:])
    assert_array_equal(brc_partial.centroids_, brc.centroids_)
