import numpy as np
import scipy.sparse as sp

from .. import as_float_array, atleast2d_or_csr, safe_asarray


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
