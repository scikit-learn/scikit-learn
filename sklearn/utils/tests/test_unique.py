import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._testing import skip_if_array_api_compat_not_configured
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


@pytest.mark.parametrize("dtype", ["float64", "U", "O", "T"])
def test_operation_cache_computes_once_and_expires(dtype, monkeypatch):
    from sklearn.utils._unique import _metadata_cache

    if dtype == "T":
        pytest.importorskip("numpy", minversion="2.0")
        dtype = np.dtypes.StringDType()
    arr = np.array(["1", "2", "1"], dtype=dtype)
    original_unique = np.unique
    calls = []

    def counted_unique(*args, **kwargs):
        calls.append(1)
        return original_unique(*args, **kwargs)

    monkeypatch.setattr(np, "unique", counted_unique)
    with _metadata_cache():
        attached = attach_unique(arr)
        first = cached_unique(attached)
        with _metadata_cache():
            assert cached_unique(attached) is first
            assert cached_unique(arr) is first
            assert attach_unique(attached) is attached
        assert len(calls) == 1
    assert arr.dtype.metadata is None
    arr[0] = "3"
    with _metadata_cache():
        assert_array_equal(cached_unique(arr), np.array(["1", "2", "3"], dtype=dtype))
    assert len(calls) == 2


def test_operation_metadata_cleanup_on_exception():
    import gc
    import weakref

    from sklearn.utils._unique import _metadata_cache

    with pytest.raises(RuntimeError, match="test"):
        with _metadata_cache():
            arr = np.array([1, 2, 1])
            reference = weakref.ref(arr)
            cached_unique(arr)
            del arr
            assert reference() is not None
            raise RuntimeError("test")
    gc.collect()
    assert reference() is None


def test_operation_metadata_does_not_reuse_slices():
    from sklearn.utils._unique import _metadata_cache

    arr = np.array([1, 2, 3])
    with _metadata_cache():
        assert_array_equal(cached_unique(arr), [1, 2, 3])
        assert_array_equal(cached_unique(arr[:1]), [1])


def test_operation_metadata_is_thread_local():
    from concurrent.futures import ThreadPoolExecutor
    from threading import Barrier

    from sklearn.utils._unique import _metadata_cache

    barrier = Barrier(2)
    arr = np.array([1, 2, 1])

    def run():
        with _metadata_cache():
            result = cached_unique(arr)
            barrier.wait(timeout=10)
            assert cached_unique(arr) is result
            return result

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(run) for _ in range(2)]
        first, second = [future.result() for future in futures]
    assert first is not second
    assert_array_equal(first, second)


@pytest.mark.parametrize("column", [False, True])
def test_stringdtype_metric_reuses_unique_after_reshape(column, monkeypatch):
    from sklearn.metrics import accuracy_score

    pytest.importorskip("numpy", minversion="2.0")
    true = np.array(["a", "b", "a"], dtype=np.dtypes.StringDType())
    pred = np.array(["b", "b", "a"], dtype=np.dtypes.StringDType())
    if column:
        true, pred = true[:, None], pred[:, None]
    original_unique = np.unique
    calls = []

    def counted_unique(*args, **kwargs):
        calls.append(1)
        return original_unique(*args, **kwargs)

    monkeypatch.setattr(np, "unique", counted_unique)
    assert accuracy_score(true, pred) == pytest.approx(2 / 3)
    assert len(calls) == 2
    true[0] = "b"
    assert accuracy_score(true, pred) == 1
    assert len(calls) == 4


@skip_if_array_api_compat_not_configured
def test_operation_metadata_array_api(monkeypatch):
    from sklearn import config_context
    from sklearn.utils._unique import _metadata_cache

    xp = pytest.importorskip("array_api_strict")
    arr = xp.asarray([1, 2, 1])
    unique_values = xp.unique_values
    calls = []

    def counted_unique(values):
        calls.append(1)
        return unique_values(values)

    monkeypatch.setattr(xp, "unique_values", counted_unique)
    with config_context(array_api_dispatch=True), _metadata_cache():
        first = cached_unique(arr)
        assert cached_unique(arr) is first
        assert len(calls) == 1
        attached = _attach_metadata(arr, source="test")
        assert attached is arr
        assert _get_metadata(attached)["source"] == "test"
    assert _get_metadata(arr) is None


def test_operation_metadata_fallback_is_not_specific_to_stringdtype(monkeypatch):
    from sklearn.utils._unique import _metadata_cache

    arr = np.array([1, 2])

    def unsupported_metadata(*args, **kwargs):
        raise TypeError("cannot attach metadata")

    with _metadata_cache():
        monkeypatch.setattr(np, "dtype", unsupported_metadata)
        assert _attach_metadata(arr, source="test") is arr
        assert _get_metadata(arr) == {"source": "test"}


@pytest.mark.parametrize("dtype", ["int64", "U", "O", "T"])
@pytest.mark.parametrize("operation", ["targets", "fit", "fit_transform", "transform"])
def test_label_validation_reuses_unique(dtype, operation, monkeypatch):
    from sklearn.preprocessing import LabelBinarizer
    from sklearn.utils.multiclass import check_classification_targets

    if dtype == "T":
        pytest.importorskip("numpy", minversion="2.0")
        dtype = np.dtypes.StringDType()
    y = np.array(["0", "1", "2"] * 10, dtype=dtype)
    estimator = LabelBinarizer()
    if operation == "transform":
        estimator.fit(y)
    original_unique = np.unique
    calls = []

    def counted_unique(values, *args, **kwargs):
        # Ignore discovery on the much smaller array of merged classes.
        if np.asarray(values).size == y.size:
            calls.append(1)
        return original_unique(values, *args, **kwargs)

    monkeypatch.setattr(np, "unique", counted_unique)
    function = (
        check_classification_targets
        if operation == "targets"
        else getattr(estimator, operation)
    )
    function(y)
    assert len(calls) == 1
    # A separate public call must recompute, even when given the same array.
    y[0] = "1"
    result = function(y)
    assert len(calls) == 2
    if operation in ("fit_transform", "transform"):
        assert_array_equal(result[0], [0, 1, 0])


@pytest.mark.parametrize("use_scope", [False, True])
def test_operation_metadata_expires_in_copied_context(use_scope):
    import gc
    import weakref
    from contextlib import nullcontext
    from contextvars import copy_context

    from sklearn.utils._unique import _metadata_cache

    with _metadata_cache():
        original = np.array([1, 2, 1])
        reference = weakref.ref(original)
        cached_unique(original)
        copied = copy_context()
    del original
    gc.collect()
    assert reference() is None

    values = np.array([1, 2, 1])

    def read():
        with _metadata_cache() if use_scope else nullcontext():
            return cached_unique(values)

    assert_array_equal(copied.run(read), [1, 2])
    values[0] = 3
    assert_array_equal(copied.run(read), [1, 2, 3])
