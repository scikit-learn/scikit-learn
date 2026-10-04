import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._unique import (
    _attach_metadata,
    _get_metadata,
    attach_unique,
    cached_unique,
)
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


@pytest.mark.parametrize(
    "numpy_string_dtype",
    [
        "U",
        "O",
        "T",
    ],
    indirect=True,
)
def test_unique_string_dtypes(numpy_string_dtype):
    arr = np.array(["b", "a", "b"], dtype=numpy_string_dtype)
    attached = attach_unique(arr)
    assert_array_equal(attached, arr)
    assert attached.dtype == arr.dtype
    # Caching is optional for dtypes that cannot carry metadata.
    assert_array_equal(cached_unique(attached), ["a", "b"])


@pytest.mark.parametrize("dtype", [np.float64, "U2", "O"])
def test_metadata_helpers_preserve_input_and_existing_metadata(dtype, monkeypatch):
    dtype = np.dtype(dtype, metadata={"source": "original"})
    arr = np.array([1, 2, 1], dtype=dtype)
    attached = attach_unique(arr)
    assert _get_metadata(arr) == {"source": "original"}
    assert _get_metadata(attached)["source"] == "original"
    unique = _get_metadata(attached)["unique"]
    assert attached.base is arr
    assert_array_equal(unique, np.unique(arr))

    def unexpected_unique(*args, **kwargs):
        pytest.fail("The attached cache should be reused")

    monkeypatch.setattr(np, "unique", unexpected_unique)
    assert attach_unique(attached) is attached
    assert cached_unique(attached) is unique


def test_string_dtype_metadata_fallback_preserves_parameters(monkeypatch):
    pytest.importorskip("numpy", minversion="2.0")
    dtype = np.dtypes.StringDType(na_object=np.nan, coerce=False)
    arr = np.array(["b", "a", "b"], dtype=dtype)
    with monkeypatch.context() as patch:

        def unexpected_unique(*args, **kwargs):
            pytest.fail("Do not compute a cache that cannot be attached")

        patch.setattr(np, "unique", unexpected_unique)
        assert attach_unique(arr) is arr
        assert _attach_metadata(arr, source="test") is arr
    assert _get_metadata(arr) is None
    assert arr.dtype is dtype
    assert arr.dtype.coerce is False
    assert arr.dtype.na_object is np.nan
    assert_array_equal(cached_unique(arr), ["a", "b"])
    # No stale cache is left behind on the unchanged input.
    arr[0] = "c"
    assert_array_equal(cached_unique(arr), ["a", "b", "c"])
