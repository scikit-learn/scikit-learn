"""Tests for input validation functions"""

from tempfile import NamedTemporaryFile
import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sp
from nose.tools import assert_raises, assert_true, assert_false, assert_equal

from sklearn.utils import (array2d, as_float_array, atleast2d_or_csr,
                           atleast2d_or_csc, check_arrays, safe_asarray)

from sklearn.random_projection import sparse_random_matrix


def test_safe_asarray():
    """Test that array dtype conversion works."""
    # Test with sparse arrays
    X = sp.csc_matrix(np.arange(4, dtype=np.float))
    Y = safe_asarray(X)
    assert_true(Y.dtype == np.float)
    # Check that no copy has been performed
    Y.data[0] = 7  # value not in original array
    assert_equal(X.data[0], Y.data[0])

    Y = safe_asarray(X, dtype=np.int)
    assert_equal(Y.data.dtype, np.int)

    # Test with dense arrays
    X = np.arange(4, dtype=np.float)
    Y = safe_asarray(X)
    assert_true(Y.dtype == np.float)
    # Check that no copy has been performed
    Y[0] = 7
    assert_equal(X[0], Y[0])

    Y = safe_asarray(X, dtype=np.int)
    assert_equal(Y.dtype, np.int)

    # Non-regression: LIL and DOK used to fail for lack of a .data attribute
    X = np.ones([2, 3])
    safe_asarray(sp.dok_matrix(X))
    safe_asarray(sp.lil_matrix(X), dtype=X.dtype)


def test_as_float_array():
    """Test function for as_float_array"""
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

    # Test the copy parameter with some matrices
    matrices = [
        np.matrix(np.arange(5)),
        sp.csc_matrix(np.arange(5)).toarray(),
        sparse_random_matrix(10, 10, density=0.10).toarray()
    ]
    for M in matrices:
        N = as_float_array(M, copy=True)
        N[0, 0] = np.nan
        assert_false(np.isnan(M).any())


def test_atleast2d_or_sparse():
    for typ in [sp.csr_matrix, sp.dok_matrix, sp.lil_matrix, sp.coo_matrix]:
        X = typ(np.arange(9, dtype=float).reshape(3, 3))

        Y = atleast2d_or_csr(X, copy=True)
        assert_true(isinstance(Y, sp.csr_matrix))
        Y.data[:] = 1
        assert_array_equal(X.toarray().ravel(), np.arange(9))

        Y = atleast2d_or_csc(X, copy=False)
        Y.data[:] = 4
        assert_true(np.all(X.data == 4)
                    if isinstance(X, sp.csc_matrix)
                    else np.all(X.toarray().ravel() == np.arange(9)))

        Y = atleast2d_or_csr(X, dtype=np.float32)
        assert_true(Y.dtype == np.float32)


def test_check_arrays_exceptions():
    """Check that invalid arguments raise appropriate exceptions"""
    assert_raises(ValueError, check_arrays, [0], [0, 1])
    assert_raises(TypeError, check_arrays, 0, [0, 1])
    assert_raises(TypeError, check_arrays, [0], 0)
    assert_raises(TypeError, check_arrays, [0, 1], [0, 1], meaning_of_life=42)
    assert_raises(ValueError, check_arrays, [0], [0], sparse_format='fake')
    assert_raises(ValueError, check_arrays, np.zeros((2, 3, 4)), [0])


def test_np_matrix():
    """Confirm that input validation code does not return np.matrix"""
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
    """Confirm that input validation code doesn't copy memory mapped arrays"""

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
    """Check that ordering is enforced correctly by validation utilities.

    We need to check each validation utility, because a 'copy' without
    'order=K' will kill the ordering.
    """
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

    X = sp.csr_matrix(X)
    X.data = X.data[::-1]
    assert_false(X.data.flags['C_CONTIGUOUS'])

    for validator in (atleast2d_or_csc, atleast2d_or_csr):
        for copy in (True, False):
            Y = validator(X, copy=copy, order='C')
            assert_true(Y.data.flags['C_CONTIGUOUS'])


def test_check_arrays():
    # check that error is raised on different length inputs
    X = [0, 1]
    Y = np.arange(3)
    assert_raises(ValueError, check_arrays, X, Y)

    # check error for sparse matrix and array
    X = sp.csc_matrix(np.arange(4))
    assert_raises(ValueError, check_arrays, X, Y)

    # check they y=None pattern
    X = [0, 1, 2]
    X_, Y_, Z_ = check_arrays(X, Y, None)
    assert_true(Z_ is None)

    # check that lists are converted
    X_, Y_ = check_arrays(X, Y)
    assert_true(isinstance(X_, np.ndarray))
    assert_true(isinstance(Y_, np.ndarray))

    # check that Y was not copied:
    assert_true(Y_ is Y)

    # check copying
    X_, Y_ = check_arrays(X, Y, copy=True)
    assert_false(Y_ is Y)

    # check forcing dtype
    X_, Y_ = check_arrays(X, Y, dtype=np.int)
    assert_equal(X_.dtype, np.int)
    assert_equal(Y_.dtype, np.int)

    X_, Y_ = check_arrays(X, Y, dtype=np.float)
    assert_equal(X_.dtype, np.float)
    assert_equal(Y_.dtype, np.float)

    # test check_ccontiguous
    Y = np.arange(6).reshape(3, 2).copy('F')
    # if we don't specify it, it is not changed
    X_, Y_ = check_arrays(X, Y)
    assert_true(Y_.flags['F_CONTIGUOUS'])
    assert_false(Y_.flags['C_CONTIGUOUS'])

    X_, Y_ = check_arrays(X, Y, check_ccontiguous=True)
    assert_true(Y_.flags['C_CONTIGUOUS'])
    assert_false(Y_.flags['F_CONTIGUOUS'])

    # check that lists are passed through if allow_lists is true
    X_, Y_ = check_arrays(X, Y, allow_lists=True)
    assert_true(isinstance(X_, list))
