"""
Tests for input validation functions
"""

from tempfile import NamedTemporaryFile
import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sp
from nose.tools import assert_raises, assert_true, assert_false

from sklearn.utils import (array2d, as_float_array, atleast2d_or_csr,
                           atleast2d_or_csc, check_arrays, safe_asarray)


def test_as_float_array():
    """
    Test function for as_float_array
    """
    X = np.ones((3, 10), dtype=np.int32)
    X = X + np.arange(10, dtype=np.int32)
    # Checks that the return type is ok
    X2 = as_float_array(X, copy=False)
    np.testing.assert_equal(X2.dtype, np.float32)
    # Another test
    X = X.astype(np.int64)
    X2 = as_float_array(X, copy=True)
    # Checking that the array wasn't overwritten
    assert_true(as_float_array(X, False) is not X)
    # Checking that the new type is ok
    np.testing.assert_equal(X2.dtype, np.float64)
    # Here, X is of the right type, it shouldn't be modified
    X = np.ones((3, 2), dtype=np.float32)
    assert_true(as_float_array(X, copy=False) is X)
    # Test that if X is fortran ordered it stays
    X = np.asfortranarray(X)
    assert_true(np.isfortran(as_float_array(X, copy=True)))


def test_check_arrays_exceptions():
    """Check that invalid arguments raise appropriate exceptions"""
    assert_raises(ValueError, check_arrays, [0], [0, 1])
    assert_raises(TypeError, check_arrays, 0, [0, 1])
    assert_raises(TypeError, check_arrays, [0], 0)
    assert_raises(TypeError, check_arrays, [0, 1], [0, 1], meaning_of_life=42)
    assert_raises(ValueError, check_arrays, [0], [0], sparse_format='fake')


def test_np_matrix():
    """
    Confirm that input validation code does not return np.matrix
    """
    X = np.arange(12).reshape(3, 4)

    assert_false(isinstance(as_float_array(X), np.matrix))
    assert_false(isinstance(as_float_array(np.matrix(X)), np.matrix))
    assert_false(isinstance(as_float_array(sp.csc_matrix(X)), np.matrix))

    assert_false(isinstance(atleast2d_or_csr(X), np.matrix))
    assert_false(isinstance(atleast2d_or_csr(np.matrix(X)), np.matrix))
    assert_false(isinstance(atleast2d_or_csr(sp.csc_matrix(X)), np.matrix))

    assert_false(isinstance(atleast2d_or_csc(X), np.matrix))
    assert_false(isinstance(atleast2d_or_csc(np.matrix(X)), np.matrix))
    assert_false(isinstance(atleast2d_or_csc(sp.csr_matrix(X)), np.matrix))

    assert_false(isinstance(safe_asarray(X), np.matrix))
    assert_false(isinstance(safe_asarray(np.matrix(X)), np.matrix))
    assert_false(isinstance(safe_asarray(sp.lil_matrix(X)), np.matrix))

    assert_true(atleast2d_or_csr(X, copy=False) is X)
    assert_false(atleast2d_or_csr(X, copy=True) is X)
    assert_true(atleast2d_or_csc(X, copy=False) is X)
    assert_false(atleast2d_or_csc(X, copy=True) is X)


def test_memmap():
    """
    Confirm that input validation code doesn't copy memory mapped arrays
    """

    asflt = lambda x: as_float_array(x, copy=False)

    with NamedTemporaryFile(prefix='sklearn-test') as tmp:
        M = np.memmap(tmp, shape=100, dtype=np.float32)
        M[:] = 0

        for f in (array2d, np.asarray, asflt, safe_asarray):
            X = f(M)
            X[:] = 1
            assert_array_equal(X.ravel(), M)
            X[:] = 0


def test_ordering():
    # Check that ordering is enforced correctly by the different
    # validation utilities
    # We need to check each validation utility, because a 'copy' without
    # 'order=K' will kill the ordering
    X = np.ones((10, 5))
    for A in X, X.T:
        for validator in (array2d, atleast2d_or_csr, atleast2d_or_csc):
            for copy in (True, False):
                B = validator(A, order='C', copy=copy)
                assert_true(B.flags['C_CONTIGUOUS'])
                B = validator(A, order='F', copy=copy)
                assert_true(B.flags['F_CONTIGUOUS'])
                if copy:
                    assert_false(A is B)
