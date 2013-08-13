"""Tests for bicluster utilities."""

import numpy as np

from scipy.sparse import csr_matrix, issparse

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_true

from sklearn.cluster.bicluster.utils import get_indicators
from sklearn.cluster.bicluster.utils import get_shape
from sklearn.cluster.bicluster.utils import get_submatrix


def test_get_indicators():
    rows = [2, 4, 5]
    columns = [0, 1, 3]
    shape = (6, 4)
    row_ind, col_ind = get_indicators(rows, columns, shape)
    assert_array_equal(row_ind, [False, False, True, False, True, True])
    assert_array_equal(col_ind, [True, True, False, True])


def test_get_shape():
    rows = [True, True, False, False]
    cols = [True, False, True, True]
    assert_equal(get_shape(rows, cols), (2, 3))


def test_get_submatrix():
    data = np.arange(20).reshape(5, 4)
    rows = [True, True, False, False, True]
    cols = [False, False, True, True]
    for X in (data, csr_matrix(data)):
        submatrix = get_submatrix(rows, cols, X)
        if issparse(submatrix):
            submatrix = submatrix.todense()
        assert_array_equal(submatrix, [[2, 3],
                                       [6, 7],
                                       [18, 19]])
        submatrix[:] = -1
        assert_true(np.all(X != -1))
