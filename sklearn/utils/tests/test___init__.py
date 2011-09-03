import numpy as np
from .. import as_float_array


def test_as_float_array():
    """
    Test function for as_float_array
    """
    X = np.ones((3, 10), dtype=np.int32)
    X = X + np.arange(10, dtype=np.int32)
    # Checks that the return type is ok
    X2 = as_float_array(X, overwrite_X=True)
    np.testing.assert_equal(X2.dtype, np.float32)
    # Another test
    X = X.astype(np.int64)
    X2 = as_float_array(X, overwrite_X=False)
    # Checking that the array wasn't overwritten
    assert as_float_array(X, False) is not X
    # Checking that the new type is ok
    np.testing.assert_equal(X2.dtype, np.float64)
    # Here, X is of the right type, it shouldn't be modified
    X = np.ones((3, 2), dtype=np.float32)
    assert as_float_array(X, overwrite_X=True) is X
