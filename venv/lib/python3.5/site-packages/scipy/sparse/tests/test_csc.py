from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_
from scipy.sparse import csr_matrix, csc_matrix


def test_csc_getrow():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsc = csc_matrix(X)

    for i in range(N):
        arr_row = X[i:i + 1, :]
        csc_row = Xcsc.getrow(i)

        assert_array_almost_equal(arr_row, csc_row.toarray())
        assert_(type(csc_row) is csr_matrix)


def test_csc_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsc = csc_matrix(X)

    for i in range(N):
        arr_col = X[:, i:i + 1]
        csc_col = Xcsc.getcol(i)

        assert_array_almost_equal(arr_col, csc_col.toarray())
        assert_(type(csc_col) is csc_matrix)

