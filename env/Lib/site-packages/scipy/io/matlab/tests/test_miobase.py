""" Testing miobase module
"""

import numpy as np

from numpy.testing import assert_equal
from pytest import raises as assert_raises

from scipy.io.matlab.miobase import matdims


def test_matdims():
    # Test matdims dimension finder
    assert_equal(matdims(np.array(1)), (1, 1))  # numpy scalar
    assert_equal(matdims(np.array([1])), (1, 1))  # 1d array, 1 element
    assert_equal(matdims(np.array([1,2])), (2, 1))  # 1d array, 2 elements
    assert_equal(matdims(np.array([[2],[3]])), (2, 1))  # 2d array, column vector
    assert_equal(matdims(np.array([[2,3]])), (1, 2))  # 2d array, row vector
    # 3d array, rowish vector
    assert_equal(matdims(np.array([[[2,3]]])), (1, 1, 2))
    assert_equal(matdims(np.array([])), (0, 0))  # empty 1d array
    assert_equal(matdims(np.array([[]])), (0, 0))  # empty 2d
    assert_equal(matdims(np.array([[[]]])), (0, 0, 0))  # empty 3d
    # Optional argument flips 1-D shape behavior.
    assert_equal(matdims(np.array([1,2]), 'row'), (1, 2))  # 1d array, 2 elements
    # The argument has to make sense though
    assert_raises(ValueError, matdims, np.array([1,2]), 'bizarre')
    # Check empty sparse matrices get their own shape
    from scipy.sparse import csr_matrix, csc_matrix
    assert_equal(matdims(csr_matrix(np.zeros((3, 3)))), (3, 3))
    assert_equal(matdims(csc_matrix(np.zeros((2, 2)))), (2, 2))
