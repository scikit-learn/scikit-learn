# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause

import numpy as np
from scipy import sparse

from sklearn.utils.graph import graph_laplacian


def test_graph_laplacian():
    for mat in (np.arange(10) * np.arange(10)[:, np.newaxis],
                np.ones((7, 7)),
                np.eye(19),
                np.vander(np.arange(4)) + np.vander(np.arange(4)).T,):
        sp_mat = sparse.csr_matrix(mat)
        for normed in (True, False):
            laplacian = graph_laplacian(mat, normed=normed)
            n_nodes = mat.shape[0]
            if not normed:
                np.testing.assert_array_almost_equal(laplacian.sum(axis=0),
                                                     np.zeros(n_nodes))
            np.testing.assert_array_almost_equal(laplacian.T, laplacian)
            np.testing.assert_array_almost_equal(
                laplacian, graph_laplacian(sp_mat, normed=normed).toarray())


def test_graph_laplacian_with_copy():
    for mat in (np.arange(10) * np.arange(10)[:, np.newaxis],
                np.ones((7, 7)),
                np.eye(19),
                np.vander(np.arange(4)) + np.vander(np.arange(4)).T,):
        sparse_mat = sparse.csr_matrix(mat)
        for normed in (True, False):
            laplacian = graph_laplacian(mat, normed=normed)
            laplacian_no_copy = graph_laplacian(mat, normed=normed, copy=False)
            if not normed:
                assert(mat is laplacian_no_copy)
            np.testing.assert_array_equal(laplacian, laplacian_no_copy)
            sparse_laplacian = graph_laplacian(sparse_mat, normed=normed)
            sparse_laplacian_no_copy = graph_laplacian(sparse_mat,
                                                       normed=normed,
                                                       copy=False)
            assert((sparse_laplacian != sparse_laplacian_no_copy).nnz == 0)
