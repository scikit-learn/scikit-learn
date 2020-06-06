from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_
from scipy.sparse import csr_matrix, csc_matrix

import pytest


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

@pytest.mark.parametrize("matrix_input, axis, expected_shape",
    [(csc_matrix([[1, 0],
                [0, 0],
                [0, 2]]),
      0, (0, 2)),
     (csc_matrix([[1, 0],
                [0, 0],
                [0, 2]]),
      1, (3, 0)),
     (csc_matrix([[1, 0],
                [0, 0],
                [0, 2]]),
      'both', (0, 0)),
     (csc_matrix([[0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 2, 3, 0, 1]]),
      0, (0, 6))])
def test_csc_empty_slices(matrix_input, axis, expected_shape):
    # see gh-11127 for related discussion
    slice_1 = matrix_input.A.shape[0] - 1
    slice_2 = slice_1
    slice_3 = slice_2 - 1

    if axis == 0:
        actual_shape_1 = matrix_input[slice_1:slice_2, :].A.shape
        actual_shape_2 = matrix_input[slice_1:slice_3, :].A.shape
    elif axis == 1:
        actual_shape_1 = matrix_input[:, slice_1:slice_2].A.shape
        actual_shape_2 = matrix_input[:, slice_1:slice_3].A.shape
    elif axis == 'both':
        actual_shape_1 = matrix_input[slice_1:slice_2, slice_1:slice_2].A.shape
        actual_shape_2 = matrix_input[slice_1:slice_3, slice_1:slice_3].A.shape

    assert actual_shape_1 == expected_shape
    assert actual_shape_1 == actual_shape_2
