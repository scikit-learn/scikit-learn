"""
Tests for the birch clusterer.
"""

import numpy as np

from sklearn.cluster.birch import Birch
from sklearn.datasets import make_blobs

from sklearn.utils.testing import assert_equal


def test_n_samples_leaves_roots():
    """Sanity check for the number of samples in leaves and roots"""
    X, y = make_blobs()
    brc = Birch()
    brc.fit(X)
    n_samples_root = sum([sc.n_ for sc in brc.root_.subclusters_])
    n_samples_leaves = sum([sc.n_ for leaf in brc.get_leaves()
                            for sc in leaf.subclusters_])
    assert_equal(n_samples_leaves, X.shape[0])
    assert_equal(n_samples_root, X.shape[0])
