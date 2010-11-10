"""
Testing for Clustering methods

"""

import numpy as np
from numpy.testing import assert_equal
from scipy import sparse
import nose

from .. import SpectralClustering


def test_spectral_clustering():
    S = np.array([[1, 5, 2, 1, 0, 0, 0],
                  [5, 1, 3, 1, 0, 0, 0],
                  [2, 3, 1, 1, 0, 0, 0],
                  [1, 1, 1, 1, 2, 1, 1],
                  [0, 0, 0, 2, 2, 3, 2],
                  [0, 0, 0, 1, 3, 1, 4],
                  [0, 0, 0, 1, 2, 4, 1],
                 ])

    for mat in (S, sparse.csr_matrix(S)):
        labels = SpectralClustering().fit(mat, k=2).labels_
        if labels[0] == 0:
            labels = 1 - labels

        assert_equal(labels, [1, 1, 1, 0, 0, 0, 0])


def test_spectral_clustering_sparse():
    # We need a large matrice, or the lobpcg solver will fallback to its
    # non-sparse and buggy mode
    raise nose.SkipTest("XFailed Test")
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

    labels = SpectralClustering().fit(S, k=2).labels_
    if labels[0] == 0:
        labels = 1 - labels

    assert np.mean(labels == [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]) > .9

