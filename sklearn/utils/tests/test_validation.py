"""
Tests for input validation functions
"""

import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sp
from tempfile import NamedTemporaryFile
from nose.tools import assert_raises

from .. import (array2d, as_float_array, atleast2d_or_csr, check_arrays,
                safe_asarray)


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
    assert as_float_array(X, False) is not X
    # Checking that the new type is ok
    np.testing.assert_equal(X2.dtype, np.float64)
    # Here, X is of the right type, it shouldn't be modified
    X = np.ones((3, 2), dtype=np.float32)
    assert as_float_array(X, copy=False) is X


def test_check_arrays_exceptions():
    """Check that invalid arguments raise appropriate exceptions"""
    assert_raises(ValueError, check_arrays, [0], [0, 1])
    assert_raises(TypeError, check_arrays, 0, [0, 1])
    assert_raises(TypeError, check_arrays, [0], 0)
    assert_raises(ValueError, check_arrays, [0, 1], [0, 1], meaning_of_life=42)
    assert_raises(ValueError, check_arrays, [0], [0], sparse_format='fake')


def test_np_matrix():
    """
    Confirm that input validation code does not return np.matrix
    """
    X = np.arange(12).reshape(3, 4)

    assert not isinstance(as_float_array(X), np.matrix)
    assert not isinstance(as_float_array(np.matrix(X)), np.matrix)
    assert not isinstance(as_float_array(sp.csc_matrix(X)), np.matrix)

    assert not isinstance(atleast2d_or_csr(X), np.matrix)
    assert not isinstance(atleast2d_or_csr(np.matrix(X)), np.matrix)
    assert not isinstance(atleast2d_or_csr(sp.csc_matrix(X)), np.matrix)

    assert not isinstance(safe_asarray(X), np.matrix)
    assert not isinstance(safe_asarray(np.matrix(X)), np.matrix)
    assert not isinstance(safe_asarray(sp.lil_matrix(X)), np.matrix)


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
