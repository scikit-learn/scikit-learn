import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._unique import attach_unique, cached_unique
from sklearn.utils.validation import check_array


def test_attach_unique_attaches_unique_to_array():
    arr = np.array([1, 2, 2, 3, 4, 4, 5])
    arr_ = attach_unique(arr)
    assert_array_equal(arr_.dtype.metadata["unique"], np.array([1, 2, 3, 4, 5]))
    assert_array_equal(arr_, arr)


def test_cached_unique_returns_cached_unique():
    my_dtype = np.dtype(np.float64, metadata={"unique": np.array([1, 2])})
    arr = np.array([1, 2, 2, 3, 4, 4, 5], dtype=my_dtype)
    assert_array_equal(cached_unique(arr), np.array([1, 2]))


def test_attach_unique_not_ndarray():
    """Test that when not np.ndarray, we don't touch the array."""
    arr = [1, 2, 2, 3, 4, 4, 5]
    arr_ = attach_unique(arr)
    assert arr_ is arr


def test_attach_unique_returns_view():
    """Test that attach_unique returns a view of the array."""
    arr = np.array([1, 2, 2, 3, 4, 4, 5])
    arr_ = attach_unique(arr)
    assert arr_.base is arr


def test_attach_unique_return_tuple():
    """Test return_tuple argument of the function."""
    arr = np.array([1, 2, 2, 3, 4, 4, 5])
    arr_tuple = attach_unique(arr, return_tuple=True)
    assert isinstance(arr_tuple, tuple)
    assert len(arr_tuple) == 1
    assert_array_equal(arr_tuple[0], arr)

    arr_single = attach_unique(arr, return_tuple=False)
    assert isinstance(arr_single, np.ndarray)
    assert_array_equal(arr_single, arr)


def test_check_array_keeps_unique():
    """Test that check_array keeps the unique metadata."""
    arr = np.array([[1, 2, 2, 3, 4, 4, 5]])
    arr_ = attach_unique(arr)
    arr_ = check_array(arr_)
    assert_array_equal(arr_.dtype.metadata["unique"], np.array([1, 2, 3, 4, 5]))
    assert_array_equal(arr_, arr)
