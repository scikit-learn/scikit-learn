import numpy as np
import pytest

from sklearn.utils._testing import assert_allclose
from sklearn.utils.arrayfuncs import any_zero_row, min_pos


def test_min_pos():
    # Check that min_pos returns a positive value and that it's consistent
    # between float and double
    X = np.random.RandomState(0).randn(100)

    min_double = min_pos(X)
    min_float = min_pos(X.astype(np.float32))

    assert_allclose(min_double, min_float)
    assert min_double >= 0


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_min_pos_no_positive(dtype):
    # Check that the return value of min_pos is the maximum representable
    # value of the input dtype when all input elements are <= 0 (#19328)
    X = np.full(100, -1.0).astype(dtype, copy=False)

    assert min_pos(X) == np.finfo(dtype).max


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_any_zero_row(dtype):
    # Check that any_zero_row returns True when a row contains a zero
    # and False otherwise

    X = np.array([[1, 2, 3], [0, 1, 2]], dtype=dtype)
    assert not any_zero_row(X)

    # Make a zero row and check return True
    X[0, :] = 0
    assert any_zero_row(X)
