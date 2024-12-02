import os
from functools import partial

import numpy
import pytest
from numpy.testing import assert_allclose

from sklearn._config import config_context
from sklearn.base import BaseEstimator
from sklearn.utils._array_api import (
    _asarray_with_order,
    _atol_for_type,
    _average,
    _convert_to_numpy,
    _count_nonzero,
    _estimator_with_converted_arrays,
    _fill_or_add_to_diagonal,
    _is_numpy_namespace,
    _isin,
    _max_precision_float_dtype,
    _nanmax,
    _nanmean,
    _nanmin,
    _NumPyAPIWrapper,
    _ravel,
    device,
    get_namespace,
    get_namespace_and_device,
    indexing_dtype,
    supported_float_dtypes,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import (
    SkipTest,
    _array_api_for_tests,
    assert_array_equal,
    skip_if_array_api_compat_not_configured,
)
from sklearn.utils.fixes import _IS_32BIT, CSR_CONTAINERS, np_version, parse_version


@pytest.mark.parametrize("X", [numpy.asarray([1, 2, 3]), [1, 2, 3]])
def test_get_namespace_ndarray_default(X):
    """Check that get_namespace returns NumPy wrapper"""
    xp_out, is_array_api_compliant = get_namespace(X)
    assert isinstance(xp_out, _NumPyAPIWrapper)
    assert not is_array_api_compliant


def test_get_namespace_ndarray_creation_device():
    """Check expected behavior with device and creation functions."""
    X = numpy.asarray([1, 2, 3])
    xp_out, _ = get_namespace(X)

    full_array = xp_out.full(10, fill_value=2.0, device="cpu")
    assert_allclose(full_array, [2.0] * 10)

    with pytest.raises(ValueError, match="Unsupported device"):
        xp_out.zeros(10, device="cuda")


@skip_if_array_api_compat_not_configured
def test_get_namespace_ndarray_with_dispatch():
    """Test get_namespace on NumPy ndarrays."""
    array_api_compat = pytest.importorskip("array_api_compat")
    if parse_version(array_api_compat.__version__) < parse_version("1.9"):
        pytest.skip(
            reason="array_api_compat was temporarily reporting NumPy as API compliant "
            "and this test would fail"
        )

    X_np = numpy.asarray([[1, 2, 3]])

    with config_context(array_api_dispatch=True):
        xp_out, is_array_api_compliant = get_namespace(X_np)
        assert is_array_api_compliant

        # In the future, NumPy should become API compliant library and we should have
        # assert xp_out is numpy
        assert xp_out is array_api_compat.numpy


@skip_if_array_api_compat_not_configured
def test_get_namespace_array_api(monkeypatch):
    """Test get_namespace for ArrayAPI arrays."""
    xp = pytest.importorskip("array_api_strict")

    X_np = numpy.asarray([[1, 2, 3]])
    X_xp = xp.asarray(X_np)
    with config_context(array_api_dispatch=True):
        xp_out, is_array_api_compliant = get_namespace(X_xp)
        assert is_array_api_compliant

        with pytest.raises(TypeError):
            xp_out, is_array_api_compliant = get_namespace(X_xp, X_np)

        def mock_getenv(key):
            if key == "SCIPY_ARRAY_API":
                return "0"

        monkeypatch.setattr("os.environ.get", mock_getenv)
        assert os.environ.get("SCIPY_ARRAY_API") != "1"
        with pytest.raises(
            RuntimeError,
            match="scipy's own support is not enabled.",
        ):
            get_namespace(X_xp)


@pytest.mark.parametrize("array_api", ["numpy", "array_api_strict"])
def test_asarray_with_order(array_api):
    """Test _asarray_with_order passes along order for NumPy arrays."""
    xp = pytest.importorskip(array_api)

    X = xp.asarray([1.2, 3.4, 5.1])
    X_new = _asarray_with_order(X, order="F", xp=xp)

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["F_CONTIGUOUS"]


@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize(
    "weights, axis, normalize, expected",
    [
        # normalize = True
        (None, None, True, 3.5),
        (None, 0, True, [2.5, 3.5, 4.5]),
        (None, 1, True, [2, 5]),
        ([True, False], 0, True, [1, 2, 3]),  # boolean weights
        ([True, True, False], 1, True, [1.5, 4.5]),  # boolean weights
        ([0.4, 0.1], 0, True, [1.6, 2.6, 3.6]),
        ([0.4, 0.2, 0.2], 1, True, [1.75, 4.75]),
        ([1, 2], 0, True, [3, 4, 5]),
        ([1, 1, 2], 1, True, [2.25, 5.25]),
        ([[1, 2, 3], [1, 2, 3]], 0, True, [2.5, 3.5, 4.5]),
        ([[1, 2, 1], [2, 2, 2]], 1, True, [2, 5]),
        # normalize = False
        (None, None, False, 21),
        (None, 0, False, [5, 7, 9]),
        (None, 1, False, [6, 15]),
        ([True, False], 0, False, [1, 2, 3]),  # boolean weights
        ([True, True, False], 1, False, [3, 9]),  # boolean weights
        ([0.4, 0.1], 0, False, [0.8, 1.3, 1.8]),
        ([0.4, 0.2, 0.2], 1, False, [1.4, 3.8]),
        ([1, 2], 0, False, [9, 12, 15]),
        ([1, 1, 2], 1, False, [9, 21]),
        ([[1, 2, 3], [1, 2, 3]], 0, False, [5, 14, 27]),
        ([[1, 2, 1], [2, 2, 2]], 1, False, [8, 30]),
    ],
)
def test_average(
    array_namespace, device_, dtype_name, weights, axis, normalize, expected
):
    xp = _array_api_for_tests(array_namespace, device_)
    array_in = numpy.asarray([[1, 2, 3], [4, 5, 6]], dtype=dtype_name)
    array_in = xp.asarray(array_in, device=device_)
    if weights is not None:
        weights = numpy.asarray(weights, dtype=dtype_name)
        weights = xp.asarray(weights, device=device_)

    with config_context(array_api_dispatch=True):
        result = _average(array_in, axis=axis, weights=weights, normalize=normalize)

    if np_version < parse_version("2.0.0") or np_version >= parse_version("2.1.0"):
        # NumPy 2.0 has a problem with the device attribute of scalar arrays:
        # https://github.com/numpy/numpy/issues/26850
        assert device(array_in) == device(result)

    result = _convert_to_numpy(result, xp)
    assert_allclose(result, expected, atol=_atol_for_type(dtype_name))


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(include_numpy_namespaces=False),
)
def test_average_raises_with_wrong_dtype(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)

    array_in = numpy.asarray([2, 0], dtype=dtype_name) + 1j * numpy.asarray(
        [4, 3], dtype=dtype_name
    )
    complex_type_name = array_in.dtype.name
    if not hasattr(xp, complex_type_name):
        # This is the case for cupy as of March 2024 for instance.
        pytest.skip(f"{array_namespace} does not support {complex_type_name}")

    array_in = xp.asarray(array_in, device=device)

    err_msg = "Complex floating point values are not supported by average."
    with (
        config_context(array_api_dispatch=True),
        pytest.raises(NotImplementedError, match=err_msg),
    ):
        _average(array_in)


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(include_numpy_namespaces=True),
)
@pytest.mark.parametrize(
    "axis, weights, error, error_msg",
    (
        (
            None,
            [1, 2],
            TypeError,
            "Axis must be specified",
        ),
        (
            0,
            [[1, 2]],
            # NumPy 2 raises ValueError, NumPy 1 raises TypeError
            (ValueError, TypeError),
            "weights",  # the message is different for NumPy 1 and 2...
        ),
        (
            0,
            [1, 2, 3, 4],
            ValueError,
            "weights",
        ),
        (0, [-1, 1], ZeroDivisionError, "Weights sum to zero, can't be normalized"),
    ),
)
def test_average_raises_with_invalid_parameters(
    array_namespace, device, dtype_name, axis, weights, error, error_msg
):
    xp = _array_api_for_tests(array_namespace, device)

    array_in = numpy.asarray([[1, 2, 3], [4, 5, 6]], dtype=dtype_name)
    array_in = xp.asarray(array_in, device=device)

    weights = numpy.asarray(weights, dtype=dtype_name)
    weights = xp.asarray(weights, device=device)

    with config_context(array_api_dispatch=True), pytest.raises(error, match=error_msg):
        _average(array_in, axis=axis, weights=weights)


def test_device_none_if_no_input():
    assert device() is None

    assert device(None, "name") is None


def test_device_inspection():
    class Device:
        def __init__(self, name):
            self.name = name

        def __eq__(self, device):
            return self.name == device.name

        def __hash__(self):
            raise TypeError("Device object is not hashable")

        def __str__(self):
            return self.name

    class Array:
        def __init__(self, device_name):
            self.device = Device(device_name)

    # Sanity check: ensure our Device mock class is non hashable, to
    # accurately account for non-hashable device objects in some array
    # libraries, because of which the `device` inspection function should'nt
    # make use of hash lookup tables (in particular, not use `set`)
    with pytest.raises(TypeError):
        hash(Array("device").device)

    # Test raise if on different devices
    err_msg = "Input arrays use different devices: cpu, mygpu"
    with pytest.raises(ValueError, match=err_msg):
        device(Array("cpu"), Array("mygpu"))

    # Test expected value is returned otherwise
    array1 = Array("device")
    array2 = Array("device")

    assert array1.device == device(array1)
    assert array1.device == device(array1, array2)
    assert array1.device == device(array1, array1, array2)


# TODO: add cupy to the list of libraries once the the following upstream issue
# has been fixed:
# https://github.com/cupy/cupy/issues/8180
@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize("library", ["numpy", "array_api_strict", "torch"])
@pytest.mark.parametrize(
    "X,reduction,expected",
    [
        ([1, 2, numpy.nan], _nanmin, 1),
        ([1, -2, -numpy.nan], _nanmin, -2),
        ([numpy.inf, numpy.inf], _nanmin, numpy.inf),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmin, axis=0),
            [1.0, 2.0, 3.0],
        ),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmin, axis=1),
            [1.0, numpy.nan, 4.0],
        ),
        ([1, 2, numpy.nan], _nanmax, 2),
        ([1, 2, numpy.nan], _nanmax, 2),
        ([-numpy.inf, -numpy.inf], _nanmax, -numpy.inf),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmax, axis=0),
            [4.0, 5.0, 6.0],
        ),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmax, axis=1),
            [3.0, numpy.nan, 6.0],
        ),
        ([1, 2, numpy.nan], _nanmean, 1.5),
        ([1, -2, -numpy.nan], _nanmean, -0.5),
        ([-numpy.inf, -numpy.inf], _nanmean, -numpy.inf),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmean, axis=0),
            [2.5, 3.5, 4.5],
        ),
        (
            [[1, 2, 3], [numpy.nan, numpy.nan, numpy.nan], [4, 5, 6.0]],
            partial(_nanmean, axis=1),
            [2.0, numpy.nan, 5.0],
        ),
    ],
)
def test_nan_reductions(library, X, reduction, expected):
    """Check NaN reductions like _nanmin and _nanmax"""
    xp = pytest.importorskip(library)

    with config_context(array_api_dispatch=True):
        result = reduction(xp.asarray(X))

    result = _convert_to_numpy(result, xp)
    assert_allclose(result, expected)


@pytest.mark.parametrize(
    "namespace, _device, _dtype", yield_namespace_device_dtype_combinations()
)
def test_ravel(namespace, _device, _dtype):
    xp = _array_api_for_tests(namespace, _device)

    array = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    array_xp = xp.asarray(array, device=_device)
    with config_context(array_api_dispatch=True):
        result = _ravel(array_xp)

    result = _convert_to_numpy(result, xp)
    expected = numpy.ravel(array, order="C")

    assert_allclose(expected, result)

    if _is_numpy_namespace(xp):
        assert numpy.asarray(result).flags["C_CONTIGUOUS"]


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize("library", ["cupy", "torch"])
def test_convert_to_numpy_gpu(library):  # pragma: nocover
    """Check convert_to_numpy for GPU backed libraries."""
    xp = pytest.importorskip(library)

    if library == "torch":
        if not xp.backends.cuda.is_built():
            pytest.skip("test requires cuda")
        X_gpu = xp.asarray([1.0, 2.0, 3.0], device="cuda")
    else:
        X_gpu = xp.asarray([1.0, 2.0, 3.0])

    X_cpu = _convert_to_numpy(X_gpu, xp=xp)
    expected_output = numpy.asarray([1.0, 2.0, 3.0])
    assert_allclose(X_cpu, expected_output)


def test_convert_to_numpy_cpu():
    """Check convert_to_numpy for PyTorch CPU arrays."""
    torch = pytest.importorskip("torch")
    X_torch = torch.asarray([1.0, 2.0, 3.0], device="cpu")

    X_cpu = _convert_to_numpy(X_torch, xp=torch)
    expected_output = numpy.asarray([1.0, 2.0, 3.0])
    assert_allclose(X_cpu, expected_output)


class SimpleEstimator(BaseEstimator):
    def fit(self, X, y=None):
        self.X_ = X
        self.n_features_ = X.shape[0]
        return self


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, converter",
    [
        ("torch", lambda array: array.cpu().numpy()),
        ("array_api_strict", lambda array: numpy.asarray(array)),
        ("cupy", lambda array: array.get()),
    ],
)
def test_convert_estimator_to_ndarray(array_namespace, converter):
    """Convert estimator attributes to ndarray."""
    xp = pytest.importorskip(array_namespace)

    X = xp.asarray([[1.3, 4.5]])
    est = SimpleEstimator().fit(X)

    new_est = _estimator_with_converted_arrays(est, converter)
    assert isinstance(new_est.X_, numpy.ndarray)


@skip_if_array_api_compat_not_configured
def test_convert_estimator_to_array_api():
    """Convert estimator attributes to ArrayAPI arrays."""
    xp = pytest.importorskip("array_api_strict")

    X_np = numpy.asarray([[1.3, 4.5]])
    est = SimpleEstimator().fit(X_np)

    new_est = _estimator_with_converted_arrays(est, lambda array: xp.asarray(array))
    assert hasattr(new_est.X_, "__array_namespace__")


def test_reshape_behavior():
    """Check reshape behavior with copy and is strict with non-tuple shape."""
    xp = _NumPyAPIWrapper()
    X = xp.asarray([[1, 2, 3], [3, 4, 5]])

    X_no_copy = xp.reshape(X, (-1,), copy=False)
    assert X_no_copy.base is X

    X_copy = xp.reshape(X, (6, 1), copy=True)
    assert X_copy.base is not X.base

    with pytest.raises(TypeError, match="shape must be a tuple"):
        xp.reshape(X, -1)


def test_get_namespace_array_api_isdtype():
    """Test isdtype implementation from _NumPyAPIWrapper."""
    xp = _NumPyAPIWrapper()

    assert xp.isdtype(xp.float32, xp.float32)
    assert xp.isdtype(xp.float32, "real floating")
    assert xp.isdtype(xp.float64, "real floating")
    assert not xp.isdtype(xp.int32, "real floating")

    for dtype in supported_float_dtypes(xp):
        assert xp.isdtype(dtype, "real floating")

    assert xp.isdtype(xp.bool, "bool")
    assert not xp.isdtype(xp.float32, "bool")

    assert xp.isdtype(xp.int16, "signed integer")
    assert not xp.isdtype(xp.uint32, "signed integer")

    assert xp.isdtype(xp.uint16, "unsigned integer")
    assert not xp.isdtype(xp.int64, "unsigned integer")

    assert xp.isdtype(xp.int64, "numeric")
    assert xp.isdtype(xp.float32, "numeric")
    assert xp.isdtype(xp.uint32, "numeric")

    assert not xp.isdtype(xp.float32, "complex floating")

    assert not xp.isdtype(xp.int8, "complex floating")
    assert xp.isdtype(xp.complex64, "complex floating")
    assert xp.isdtype(xp.complex128, "complex floating")

    with pytest.raises(ValueError, match="Unrecognized data type"):
        assert xp.isdtype(xp.int16, "unknown")


@pytest.mark.parametrize(
    "namespace, _device, _dtype", yield_namespace_device_dtype_combinations()
)
def test_indexing_dtype(namespace, _device, _dtype):
    xp = _array_api_for_tests(namespace, _device)

    if _IS_32BIT:
        assert indexing_dtype(xp) == xp.int32
    else:
        assert indexing_dtype(xp) == xp.int64


@pytest.mark.parametrize(
    "namespace, _device, _dtype", yield_namespace_device_dtype_combinations()
)
def test_max_precision_float_dtype(namespace, _device, _dtype):
    xp = _array_api_for_tests(namespace, _device)
    expected_dtype = xp.float32 if _device == "mps" else xp.float64
    assert _max_precision_float_dtype(xp, _device) == expected_dtype


@pytest.mark.parametrize(
    "array_namespace, device, _", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize("invert", [True, False])
@pytest.mark.parametrize("assume_unique", [True, False])
@pytest.mark.parametrize("element_size", [6, 10, 14])
@pytest.mark.parametrize("int_dtype", ["int16", "int32", "int64", "uint8"])
def test_isin(
    array_namespace, device, _, invert, assume_unique, element_size, int_dtype
):
    xp = _array_api_for_tests(array_namespace, device)
    r = element_size // 2
    element = 2 * numpy.arange(element_size).reshape((r, 2)).astype(int_dtype)
    test_elements = numpy.array(numpy.arange(14), dtype=int_dtype)
    element_xp = xp.asarray(element, device=device)
    test_elements_xp = xp.asarray(test_elements, device=device)
    expected = numpy.isin(
        element=element,
        test_elements=test_elements,
        assume_unique=assume_unique,
        invert=invert,
    )
    with config_context(array_api_dispatch=True):
        result = _isin(
            element=element_xp,
            test_elements=test_elements_xp,
            xp=xp,
            assume_unique=assume_unique,
            invert=invert,
        )

    assert_array_equal(_convert_to_numpy(result, xp=xp), expected)


def test_get_namespace_and_device():
    # Use torch as a library with custom Device objects:
    torch = pytest.importorskip("torch")
    xp_torch = pytest.importorskip("array_api_compat.torch")
    some_torch_tensor = torch.arange(3, device="cpu")
    some_numpy_array = numpy.arange(3)

    # When dispatch is disabled, get_namespace_and_device should return the
    # default NumPy wrapper namespace and "cpu" device. Our code will handle such
    # inputs via the usual __array__ interface without attempting to dispatch
    # via the array API.
    namespace, is_array_api, device = get_namespace_and_device(some_torch_tensor)
    assert namespace is get_namespace(some_numpy_array)[0]
    assert not is_array_api
    assert device.type == "cpu"

    # Otherwise, expose the torch namespace and device via array API compat
    # wrapper.
    with config_context(array_api_dispatch=True):
        namespace, is_array_api, device = get_namespace_and_device(some_torch_tensor)
        assert namespace is xp_torch
        assert is_array_api
        assert device == some_torch_tensor.device


@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
@pytest.mark.parametrize("axis", [0, 1, None, -1, -2])
@pytest.mark.parametrize("sample_weight_type", [None, "int", "float"])
def test_count_nonzero(
    array_namespace, device_, dtype_name, csr_container, axis, sample_weight_type
):
    from sklearn.utils.sparsefuncs import count_nonzero as sparse_count_nonzero

    xp = _array_api_for_tests(array_namespace, device_)
    array = numpy.array([[0, 3, 0], [2, -1, 0], [0, 0, 0], [9, 8, 7], [4, 0, 5]])
    if sample_weight_type == "int":
        sample_weight = numpy.asarray([1, 2, 2, 3, 1])
    elif sample_weight_type == "float":
        sample_weight = numpy.asarray([0.5, 1.5, 0.8, 3.2, 2.4], dtype=dtype_name)
    else:
        sample_weight = None
    expected = sparse_count_nonzero(
        csr_container(array), axis=axis, sample_weight=sample_weight
    )
    array_xp = xp.asarray(array, device=device_)

    with config_context(array_api_dispatch=True):
        result = _count_nonzero(
            array_xp, axis=axis, sample_weight=sample_weight, xp=xp, device=device_
        )

    assert_allclose(_convert_to_numpy(result, xp=xp), expected)

    if np_version < parse_version("2.0.0") or np_version >= parse_version("2.1.0"):
        # NumPy 2.0 has a problem with the device attribute of scalar arrays:
        # https://github.com/numpy/numpy/issues/26850
        assert device(array_xp) == device(result)


@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize("wrap", [True, False])
def test_fill_or_add_to_diagonal(array_namespace, device_, dtype_name, wrap):
    xp = _array_api_for_tests(array_namespace, device_)
    array_np = numpy.zeros((5, 4), dtype=numpy.int64)
    array_xp = xp.asarray(array_np)
    _fill_or_add_to_diagonal(array_xp, value=1, xp=xp, add_value=False, wrap=wrap)
    numpy.fill_diagonal(array_np, val=1, wrap=wrap)
    assert_array_equal(_convert_to_numpy(array_xp, xp=xp), array_np)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
@pytest.mark.parametrize("dispatch", [True, False])
def test_sparse_device(csr_container, dispatch):
    a, b = csr_container(numpy.array([[1]])), csr_container(numpy.array([[2]]))
    try:
        with config_context(array_api_dispatch=dispatch):
            assert device(a, b) is None
            assert device(a, numpy.array([1])) == "cpu"
            assert get_namespace_and_device(a, b)[2] is None
            assert get_namespace_and_device(a, numpy.array([1]))[2] == "cpu"
    except ImportError:
        raise SkipTest("array_api_compat is not installed")
