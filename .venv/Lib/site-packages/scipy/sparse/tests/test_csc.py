import operator
import sys
import platform

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, csr_array, csc_array

import pytest


LINUX_INTEL = (sys.platform == 'linux') and (platform.machine() == 'x86_64')



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
    slice_1 = matrix_input.toarray().shape[0] - 1
    slice_2 = slice_1
    slice_3 = slice_2 - 1

    if axis == 0:
        actual_shape_1 = matrix_input[slice_1:slice_2, :].toarray().shape
        actual_shape_2 = matrix_input[slice_1:slice_3, :].toarray().shape
    elif axis == 1:
        actual_shape_1 = matrix_input[:, slice_1:slice_2].toarray().shape
        actual_shape_2 = matrix_input[:, slice_1:slice_3].toarray().shape
    elif axis == 'both':
        actual_shape_1 = matrix_input[slice_1:slice_2, slice_1:slice_2].toarray().shape
        actual_shape_2 = matrix_input[slice_1:slice_3, slice_1:slice_3].toarray().shape

    assert actual_shape_1 == expected_shape
    assert actual_shape_1 == actual_shape_2


@pytest.mark.parametrize('ax', (-2, -1, 0, 1, None))
def test_argmax_overflow(ax):
    # See gh-13646: Windows integer overflow for large sparse matrices.
    dim = (100000, 100000)
    A = lil_matrix(dim)
    A[-2, -2] = 42
    A[-3, -3] = 0.1234
    A = csc_matrix(A)
    idx = A.argmax(axis=ax)

    if ax is None:
        # idx is a single flattened index
        # that we need to convert to a 2d index pair;
        # can't do this with np.unravel_index because
        # the dimensions are too large
        ii = idx % dim[0]
        jj = idx // dim[0]
    else:
        # idx is an array of size of A.shape[ax];
        # check the max index to make sure no overflows
        # we encountered
        assert np.count_nonzero(idx) == A.nnz
        ii, jj = np.max(idx), np.argmax(idx)

    assert A[ii, jj] == A[-2, -2]


@pytest.mark.skipif(not LINUX_INTEL, reason="avoid variations due to OS, see gh-23826")
@pytest.mark.timeout(2)  # only slow when broken (when spurious conversion occurs)
@pytest.mark.parametrize("op", (operator.ne, operator.lt, operator.gt,
                                operator.add, operator.sub, operator.mul,
                                # truediv, eq, ge, le make dense output so not tested
                                lambda x, y: x.minimum(y), lambda x, y: x.maximum(y)))
def test_compressed_rc_conversion_mixup(op):
    # see gh-23826 for related discussion
    num_minor_axis = np.iinfo(np.uint32).max + 1
    minor_axis_index = np.array([num_minor_axis - 1])
    major_axis_index = np.array([10])
    row_cols = (minor_axis_index, major_axis_index)
    col_rows = (major_axis_index, minor_axis_index)

    X = csc_array((np.array([10]), row_cols), shape=(num_minor_axis, 20))
    X_2 = X.copy()
    # causes timeout error upon large memory alloc only if conversion to CSR occurs
    op(X_2, X)

    Z = csr_array((np.array([10]), col_rows), shape=(20, num_minor_axis))
    Z_2 = Z.copy()
    # causes timeout error upon large memory alloc only if conversion to CSC occurs
    op(Z_2, Z)
