import numpy as np
from .. import as_float_array


def test_as_float_array():
    """
    Test function for as_float_array
    """
    # Checks that the return type is ok
    X = np.ones((3, 10), dtype=np.int32)
    X = X + np.arange(10, dtype=np.int32)
    X2 = as_float_array(X, overwrite_X=True)
    np.testing.assert_equal(X2.dtype, np.float32)
    # Checking that the array was actually overwritten
    np.testing.assert_equal(X.dtype, np.float32)
    # Another test
    X = X.astype(np.int64)
    X2 = as_float_array(X, overwrite_X=False)
    # Checking that the array wasn't overwritten
    np.testing.assert_equal(X.dtype, np.int64)
    # Checking that the new type is ok
    np.testing.assert_equal(X2.dtype, np.float64)
