# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
from scipy import sparse

from ..graph import graph_laplacian


def test_graph_laplacian():
    for mat in (np.arange(10) * np.arange(10)[:, np.newaxis],
                np.ones((7, 7)),
                np.eye(19),
                np.vander(np.arange(4)) + np.vander(np.arange(4)).T,
               ):
        sp_mat = sparse.csr_matrix(mat)
        for normed in (True, False):
            laplacian = graph_laplacian(mat, normed=normed)
            n_nodes = mat.shape[0]
            if not normed:
                np.testing.assert_array_almost_equal(laplacian.sum(axis=0),
                                            np.zeros(n_nodes))
            np.testing.assert_array_almost_equal(laplacian.T,
                                        laplacian)
            np.testing.assert_array_almost_equal(laplacian,
                            graph_laplacian(sp_mat, normed=normed).todense())
