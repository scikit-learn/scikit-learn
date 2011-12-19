"""Testing for Spectral Clustering methods"""

from cPickle import dumps, loads
import nose

import numpy as np
from numpy.testing import assert_equal
from scipy import sparse

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
        model = SpectralClustering(random_state=0, k=2).fit(mat)
        labels = model.labels_
        if labels[0] == 0:
            labels = 1 - labels

        assert_equal(labels, [1, 1, 1, 0, 0, 0, 0])

        model_copy = loads(dumps(model))
        assert_equal(model_copy.k, model.k)
        assert_equal(model_copy.mode, model.mode)
        assert_equal(model_copy.random_state.get_state(),
                     model.random_state.get_state())
        assert_equal(model_copy.labels_, model.labels_)


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

    labels = SpectralClustering(random_state=0).fit(S, k=2).labels_
    if labels[0] == 0:
        labels = 1 - labels

    assert np.mean(labels == [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]) > .9
