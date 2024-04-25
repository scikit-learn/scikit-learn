import re
from functools import partial

import numpy
import pytest
from numpy.testing import assert_allclose

from sklearn._config import config_context
from sklearn.base import BaseEstimator
from sklearn.utils._array_api import (
    _ArrayAPIWrapper,
    _asarray_with_order,
    _atol_for_type,
    _average,
    _convert_to_numpy,
    _estimator_with_converted_arrays,
    _is_numpy_namespace,
    _nanmax,
    _nanmin,
    _NumPyAPIWrapper,
    _ravel,
    device,
    get_namespace,
    indexing_dtype,
    supported_float_dtypes,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import (
    _array_api_for_tests,
    skip_if_array_api_compat_not_configured,
)
from sklearn.utils.fixes import _IS_32BIT


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

    X_np = numpy.asarray([[1, 2, 3]])

    with config_context(array_api_dispatch=True):
        xp_out, is_array_api_compliant = get_namespace(X_np)
        assert is_array_api_compliant
        assert xp_out is array_api_compat.numpy


@skip_if_array_api_compat_not_configured
def test_get_namespace_array_api():
    """Test get_namespace for ArrayAPI arrays."""
    xp = pytest.importorskip("array_api_strict")

    X_np = numpy.asarray([[1, 2, 3]])
    X_xp = xp.asarray(X_np)
    with config_context(array_api_dispatch=True):
        xp_out, is_array_api_compliant = get_namespace(X_xp)
        assert is_array_api_compliant

        with pytest.raises(TypeError):
            xp_out, is_array_api_compliant = get_namespace(X_xp, X_np)


class _AdjustableNameAPITestWrapper(_ArrayAPIWrapper):
    """API wrapper that has an adjustable name. Used for testing."""

    def __init__(self, array_namespace, name):
        super().__init__(array_namespace=array_namespace)
        self.__name__ = name


def test_array_api_wrapper_astype():
    """Test _ArrayAPIWrapper for ArrayAPIs that is not NumPy."""
    array_api_strict = pytest.importorskip("array_api_strict")
    xp_ = _AdjustableNameAPITestWrapper(array_api_strict, "array_api_strict")
    xp = _ArrayAPIWrapper(xp_)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp.astype(X, xp.float32)
    assert X_converted.dtype == xp.float32

    X_converted = xp.asarray(X, dtype=xp.float32)
    assert X_converted.dtype == xp.float32


@pytest.mark.parametrize("array_api", ["numpy", "array_api_strict"])
def test_asarray_with_order(array_api):
    """Test _asarray_with_order passes along order for NumPy arrays."""
    xp = pytest.importorskip(array_api)

    X = xp.asarray([1.2, 3.4, 5.1])
    X_new = _asarray_with_order(X, order="F", xp=xp)

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["F_CONTIGUOUS"]


def test_asarray_with_order_ignored():
    """Test _asarray_with_order ignores order for Generic ArrayAPI."""
    xp = pytest.importorskip("array_api_strict")
    xp_ = _AdjustableNameAPITestWrapper(xp, "array_api_strict")

    X = numpy.asarray([[1.2, 3.4, 5.1], [3.4, 5.5, 1.2]], order="C")
    X = xp_.asarray(X)

    X_new = _asarray_with_order(X, order="F", xp=xp_)

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["C_CONTIGUOUS"]
    assert not X_new_np.flags["F_CONTIGUOUS"]


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
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
    array_namespace, device, dtype_name, weights, axis, normalize, expected
):
    xp = _array_api_for_tests(array_namespace, device)
    array_in = numpy.asarray([[1, 2, 3], [4, 5, 6]], dtype=dtype_name)
    array_in = xp.asarray(array_in, device=device)
    if weights is not None:
        weights = numpy.asarray(weights, dtype=dtype_name)
        weights = xp.asarray(weights, device=device)

    with config_context(array_api_dispatch=True):
        result = _average(array_in, axis=axis, weights=weights, normalize=normalize)

    assert getattr(array_in, "device", None) == getattr(result, "device", None)

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
            TypeError,
            "1D weights expected",
        ),
        (
            0,
            [1, 2, 3, 4],
            ValueError,
            "Length of weights",
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


def test_device_raises_if_no_input():
    err_msg = re.escape(
        "At least one input array expected after filtering with remove_none=True, "
        "remove_types=[str]. Got none. Original types: []."
    )
    with pytest.raises(ValueError, match=err_msg):
        device()

    err_msg = re.escape(
        "At least one input array expected after filtering with remove_none=True, "
        "remove_types=[str]. Got none. Original types: [NoneType, str]."
    )
    with pytest.raises(ValueError, match=err_msg):
        device(None, "name")


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


# TODO: add cupy and cupy.array_api to the list of libraries once the
# the following upstream issue has been fixed:
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
@pytest.mark.parametrize("library", ["cupy", "torch", "cupy.array_api"])
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
        ("cupy.array_api", lambda array: array._array.get()),
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


@pytest.mark.parametrize("wrapper", [_ArrayAPIWrapper, _NumPyAPIWrapper])
def test_get_namespace_array_api_isdtype(wrapper):
    """Test isdtype implementation from _ArrayAPIWrapper and _NumPyAPIWrapper."""

    if wrapper == _ArrayAPIWrapper:
        xp_ = pytest.importorskip("array_api_strict")
        xp = _ArrayAPIWrapper(xp_)
    else:
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

    if wrapper == _NumPyAPIWrapper:
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
