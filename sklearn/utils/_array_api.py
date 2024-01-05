"""Tools to support array_api."""

import itertools
import math
from functools import lru_cache, wraps

import numpy
import scipy.special as special

from .._config import get_config
from .fixes import parse_version

_NUMPY_NAMESPACE_NAMES = {"numpy", "array_api_compat.numpy", "numpy.array_api"}


def yield_namespace_device_dtype_combinations(include_numpy_namespaces=True):
    """Yield supported namespace, device, dtype tuples for testing.

    Use this to test that an estimator works with all combinations.

    Parameters
    ----------
    include_numpy_namespaces : True
        If True, also yield numpy namespaces.

    devices : list
        If not None, returns only combinations for which the device is in
        the list.

    Returns
    -------
    array_namespace : str
        The name of the Array API namespace.

    device : str
        The name of the device on which to allocate the arrays. Can be None to
        indicate that the default value should be used.

    dtype_name : str
        The name of the data type to use for arrays. Can be None to indicate
        that the default value should be used.
    """
    for array_namespace in [
        # The following is used to test the array_api_compat wrapper when
        # array_api_dispatch is enabled: in particular, the arrays used in the
        # tests are regular numpy arrays without any "device" attribute.
        "numpy",
        # Stricter NumPy-based Array API implementation. The
        # numpy.array_api.Array instances always a dummy "device" attribute.
        "numpy.array_api",
        "cupy",
        "cupy.array_api",
        "torch",
    ]:
        if not include_numpy_namespaces and array_namespace in _NUMPY_NAMESPACE_NAMES:
            continue
        if array_namespace == "torch":
            for device, dtype in itertools.product(
                ("cpu", "cuda"), ("float64", "float32")
            ):
                yield array_namespace, device, dtype
            yield array_namespace, "mps", "float32"
        else:
            yield array_namespace, None, None


def _check_array_api_dispatch(array_api_dispatch):
    """Check that array_api_compat is installed and NumPy version is compatible.

    array_api_compat follows NEP29, which has a higher minimum NumPy version than
    scikit-learn.
    """
    if array_api_dispatch:
        try:
            import array_api_compat  # noqa
        except ImportError:
            raise ImportError(
                "array_api_compat is required to dispatch arrays using the API"
                " specification"
            )

        numpy_version = parse_version(numpy.__version__)
        min_numpy_version = "1.21"
        if numpy_version < parse_version(min_numpy_version):
            raise ImportError(
                f"NumPy must be {min_numpy_version} or newer to dispatch array using"
                " the API specification"
            )


def _single_array_device(array):
    if isinstance(array, (numpy.ndarray, numpy.generic)) or not hasattr(
        array, "device"
    ):
        return "cpu"
    else:
        return array.device


def device(*array_list):
    """Hardware device the array data resides on.

    Parameters
    ----------
    *array_list : arrays
        List of array instances from NumPy or an array API compatible library.

    Returns
    -------
    out : device
        `device` object (see the "Device Support" section of the array API spec).
    """
    if not array_list:
        raise ValueError("At least one input array expected, got none.")

    device_ = _single_array_device(array_list[0])

    for array in array_list[1:]:
        device_other = _single_array_device(array)
        if device != device_other:
            raise ValueError(
                f"Input arrays use different devices: {str(device_)}, "
                f"{str(device_other)}"
            )

    return device_


def size(x):
    """Return the total number of elements of x.

    Parameters
    ----------
    x : array
        Array instance from NumPy or an array API compatible library.

    Returns
    -------
    out : int
        Total number of elements.
    """
    return math.prod(x.shape)


def _is_numpy_namespace(xp):
    """Return True if xp is backed by NumPy."""
    return xp.__name__ in _NUMPY_NAMESPACE_NAMES


def _union1d(a, b, xp):
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.union1d(a, b))
    assert a.ndim == b.ndim == 1
    return xp.unique_values(xp.concat([xp.unique_values(a), xp.unique_values(b)]))


def isdtype(dtype, kind, *, xp):
    """Returns a boolean indicating whether a provided dtype is of type "kind".

    Included in the v2022.12 of the Array API spec.
    https://data-apis.org/array-api/latest/API_specification/generated/array_api.isdtype.html
    """
    if isinstance(kind, tuple):
        return any(_isdtype_single(dtype, k, xp=xp) for k in kind)
    else:
        return _isdtype_single(dtype, kind, xp=xp)


def _isdtype_single(dtype, kind, *, xp):
    if isinstance(kind, str):
        if kind == "bool":
            return dtype == xp.bool
        elif kind == "signed integer":
            return any(
                hasattr(xp, dtype_name) and (dtype == getattr(xp, dtype_name))
                for dtype_name in ["int8", "int16", "int32", "int64"]
            )
        elif kind == "unsigned integer":
            return any(
                hasattr(xp, dtype_name) and (dtype == getattr(xp, dtype_name))
                for dtype_name in ["uint8", "uint16", "uint32", "uint64"]
            )
        elif kind == "integral":
            return any(
                _isdtype_single(dtype, k, xp=xp)
                for k in ("signed integer", "unsigned integer")
            )
        elif kind == "real floating":
            return dtype in supported_float_dtypes(xp)
        elif kind == "complex floating":
            # Some name spaces do not have complex, such as cupy.array_api
            # and numpy.array_api
            complex_dtypes = set()
            if hasattr(xp, "complex64"):
                complex_dtypes.add(xp.complex64)
            if hasattr(xp, "complex128"):
                complex_dtypes.add(xp.complex128)
            return dtype in complex_dtypes
        elif kind == "numeric":
            return any(
                _isdtype_single(dtype, k, xp=xp)
                for k in ("integral", "real floating", "complex floating")
            )
        else:
            raise ValueError(f"Unrecognized data type kind: {kind!r}")
    else:
        return dtype == kind


@lru_cache
def _supports_dtype(xp, device, dtype):
    """Check if a given namespace/device/dtype combination is supported.

    Note that some namespaces expose dtypes that can cause a failure at runtime
    when trying to allocate an array with a specific device/dtype combination.

    This is the case for the  Pytorch / mps / float64 combination:
    at the time of writing, only float16/float32 arrays can be allocated on this
    type of device.  Otherwise a `TypeError` would be raised.

    This helper function can be refactored once an expressive enough inspection
    API has been specified as part of the standard and implemented in the main
    libraries:
    https://github.com/data-apis/array-api/issues/640
    """
    if not hasattr(xp, dtype):
        return False

    dtype = getattr(xp, dtype)

    try:
        array = xp.ones((1,), device=device, dtype=dtype)
        array += array
        float(array[0])
    except Exception:
        return False

    return True


@lru_cache
def supported_float_dtypes(xp, device=None):
    """Supported floating point types for the namespace

    Note: float16 is not officially part of the Array API spec at the
    time of writing but scikit-learn estimators and functions can choose
    to accept it when xp.float16 is defined.

    https://data-apis.org/array-api/latest/API_specification/data_types.html
    """
    return tuple(
        getattr(xp, dtype)
        for dtype in ["float64", "float32", "float16"]
        if _supports_dtype(xp, device, dtype)
    )


class _ArrayAPIWrapper:
    """sklearn specific Array API compatibility wrapper

    This wrapper makes it possible for scikit-learn maintainers to
    deal with discrepancies between different implementations of the
    Python Array API standard and its evolution over time.

    The Python Array API standard specification:
    https://data-apis.org/array-api/latest/

    Documentation of the NumPy implementation:
    https://numpy.org/neps/nep-0047-array-api-standard.html
    """

    def __init__(self, array_namespace):
        self._namespace = array_namespace

    def __getattr__(self, name):
        return getattr(self._namespace, name)

    def __eq__(self, other):
        return self._namespace == other._namespace

    def __hash__(self):
        return hash((self._namespace, "_ArrayAPIWrapper"))

    def isdtype(self, dtype, kind):
        return isdtype(dtype, kind, xp=self._namespace)


def _check_device_cpu(device):  # noqa
    if device not in {"cpu", None}:
        raise ValueError(f"Unsupported device for NumPy: {device!r}")


def _accept_device_cpu(func):
    @wraps(func)
    def wrapped_func(*args, **kwargs):
        _check_device_cpu(kwargs.pop("device", None))
        return func(*args, **kwargs)

    return wrapped_func


class _NumPyAPIWrapper:
    """Array API compat wrapper for any numpy version

    NumPy < 1.22 does not expose the numpy.array_api namespace. This
    wrapper makes it possible to write code that uses the standard
    Array API while working with any version of NumPy supported by
    scikit-learn.

    See the `get_namespace()` public function for more details.
    """

    # Creation functions in spec:
    # https://data-apis.org/array-api/latest/API_specification/creation_functions.html
    _CREATION_FUNCS = {
        "arange",
        "empty",
        "empty_like",
        "eye",
        "full",
        "full_like",
        "linspace",
        "ones",
        "ones_like",
        "zeros",
        "zeros_like",
    }
    # Data types in spec
    # https://data-apis.org/array-api/latest/API_specification/data_types.html
    _DTYPES = {
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        # XXX: float16 is not part of the Array API spec but exposed by
        # some namespaces.
        "float16",
        "float32",
        "float64",
        "complex64",
        "complex128",
    }

    def __getattr__(self, name):
        attr = getattr(numpy, name)

        # Support device kwargs and make sure they are on the CPU
        if name in self._CREATION_FUNCS:
            return _accept_device_cpu(attr)

        # Convert to dtype objects
        if name in self._DTYPES:
            return numpy.dtype(attr)
        return attr

    @property
    def bool(self):
        return numpy.bool_

    def astype(self, x, dtype, *, copy=True, casting="unsafe"):
        # astype is not defined in the top level NumPy namespace
        return x.astype(dtype, copy=copy, casting=casting)

    def asarray(self, x, *, dtype=None, device=None, copy=None):  # noqa
        _check_device_cpu(device)
        # Support copy in NumPy namespace
        if copy is True:
            return numpy.array(x, copy=True, dtype=dtype)
        else:
            return numpy.asarray(x, dtype=dtype)

    def unique_inverse(self, x):
        return numpy.unique(x, return_inverse=True)

    def unique_counts(self, x):
        return numpy.unique(x, return_counts=True)

    def unique_values(self, x):
        return numpy.unique(x)

    def concat(self, arrays, *, axis=None):
        return numpy.concatenate(arrays, axis=axis)

    def reshape(self, x, shape, *, copy=None):
        """Gives a new shape to an array without changing its data.

        The Array API specification requires shape to be a tuple.
        https://data-apis.org/array-api/latest/API_specification/generated/array_api.reshape.html
        """
        if not isinstance(shape, tuple):
            raise TypeError(
                f"shape must be a tuple, got {shape!r} of type {type(shape)}"
            )

        if copy is True:
            x = x.copy()
        return numpy.reshape(x, shape)

    def isdtype(self, dtype, kind):
        return isdtype(dtype, kind, xp=self)


_NUMPY_API_WRAPPER_INSTANCE = _NumPyAPIWrapper()


def get_namespace(*arrays):
    """Get namespace of arrays.

    Introspect `arrays` arguments and return their common Array API
    compatible namespace object, if any. NumPy 1.22 and later can
    construct such containers using the `numpy.array_api` namespace
    for instance.

    See: https://numpy.org/neps/nep-0047-array-api-standard.html

    If `arrays` are regular numpy arrays, an instance of the
    `_NumPyAPIWrapper` compatibility wrapper is returned instead.

    Namespace support is not enabled by default. To enabled it
    call:

      sklearn.set_config(array_api_dispatch=True)

    or:

      with sklearn.config_context(array_api_dispatch=True):
          # your code here

    Otherwise an instance of the `_NumPyAPIWrapper`
    compatibility wrapper is always returned irrespective of
    the fact that arrays implement the `__array_namespace__`
    protocol or not.

    Parameters
    ----------
    *arrays : array objects
        Array objects.

    Returns
    -------
    namespace : module
        Namespace shared by array objects. If any of the `arrays` are not arrays,
        the namespace defaults to NumPy.

    is_array_api_compliant : bool
        True if the arrays are containers that implement the Array API spec.
        Always False when array_api_dispatch=False.
    """
    array_api_dispatch = get_config()["array_api_dispatch"]
    if not array_api_dispatch:
        return _NUMPY_API_WRAPPER_INSTANCE, False

    _check_array_api_dispatch(array_api_dispatch)

    # array-api-compat is a required dependency of scikit-learn only when
    # configuring `array_api_dispatch=True`. Its import should therefore be
    # protected by _check_array_api_dispatch to display an informative error
    # message in case it is missing.
    import array_api_compat

    namespace, is_array_api_compliant = array_api_compat.get_namespace(*arrays), True

    # These namespaces need additional wrapping to smooth out small differences
    # between implementations
    if namespace.__name__ in {
        "numpy.array_api",
        "cupy.array_api",
    }:
        namespace = _ArrayAPIWrapper(namespace)

    return namespace, is_array_api_compliant


def _expit(X):
    xp, _ = get_namespace(X)
    if _is_numpy_namespace(xp):
        return xp.asarray(special.expit(numpy.asarray(X)))

    return 1.0 / (1.0 + xp.exp(-X))


def _add_to_diagonal(array, value, xp):
    # Workaround for the lack of support for xp.reshape(a, shape, copy=False) in
    # numpy.array_api: https://github.com/numpy/numpy/issues/23410
    value = xp.asarray(value, dtype=array.dtype)
    if _is_numpy_namespace(xp):
        array_np = numpy.asarray(array)
        array_np.flat[:: array.shape[0] + 1] += value
        return xp.asarray(array_np)
    elif value.ndim == 1:
        for i in range(array.shape[0]):
            array[i, i] += value[i]
    else:
        # scalar value
        for i in range(array.shape[0]):
            array[i, i] += value


def _average(a, axis=None, weights=None, normalize=True, xp=None):
    """Partial port of np.average to support the Array API."""
    input_arrays = [a]
    if weights is not None:
        input_arrays.append(weights)

    if xp is None:
        xp, _ = get_namespace(*input_arrays)

    device_ = device(*input_arrays)

    if _is_numpy_namespace(xp) and normalize:
        return xp.asarray(numpy.average(a, axis=axis, weights=weights))

    a = xp.asarray(a, device=device_)
    if weights is not None:
        weights = xp.asarray(weights, device=device_)

    output_dtype = None
    output_dtype_name = None

    if xp.isdtype(a.dtype, "bool"):
        a = xp.astype(a, xp.int32)
    if weights is not None and xp.isdtype(weights.dtype, "bool"):
        weights = xp.astype(weights, xp.int32)

    if any(
        (
            (input_array is not None)
            and (
                (not xp.isdtype(input_array.dtype, "numeric"))
                or xp.isdtype(input_array.dtype, "complex floating")
            )
        )
        for input_array in [a, weights]
    ):
        raise ValueError("Expecting only boolean, integral or real floating values.")

    if weights is None and xp.isdtype(a.dtype, "integral"):
        output_dtype_name = "float64"
    elif weights is None:
        output_dtype = a.dtype
    elif xp.isdtype(a.dtype, "real floating") and xp.isdtype(
        weights.dtype, "real floating"
    ):
        output_dtype = (
            a.dtype
            if (xp.finfo(a.dtype).bits >= xp.finfo(a.dtype).bits)
            else weights.dtype
        )
    else:
        output_dtype_name = "float64"

    cast_to_float64 = (output_dtype_name == "float64") or (
        xp.finfo(output_dtype).bits == 64
    )

    if cast_to_float64 and not _supports_dtype(xp, device_, "float64"):
        a = _convert_to_numpy(a, xp)
        if weights is not None:
            weights = _convert_to_numpy(weights, xp)
        return xp.asarray(
            _average(a, axis=axis, normalize=normalize, weights=weights),
            dtype=xp.float32,
            device=device_,
        )

    if output_dtype is None:
        output_dtype = getattr(xp, output_dtype_name)

    a = xp.astype(a, output_dtype)

    if weights is None:
        return (xp.mean if normalize else xp.sum)(a, axis=axis)

    weights = xp.astype(weights, output_dtype)

    if a.shape != weights.shape:
        if axis is None:
            raise TypeError(
                "Axis must be specified when shapes of a and weights differ."
            )
        if weights.ndim != 1:
            raise TypeError("1D weights expected when shapes of a and weights differ.")

        if size(weights) != a.shape[axis]:
            raise ValueError("Length of weights not compatible with specified axis.")

        # If weights are 1D, add singleton dimensions for broadcasting
        shape = [1] * a.ndim
        shape[axis] = a.shape[axis]
        weights = xp.reshape(weights, shape)

    sum_ = xp.sum(xp.multiply(a, weights), axis=axis)

    if not normalize:
        return sum_

    scale = xp.sum(weights, axis=axis)
    if xp.any(scale == 0.0):
        raise ZeroDivisionError("Weights sum to zero, can't be normalized")

    return sum_ / scale


def _nanmin(X, axis=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _ = get_namespace(X)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmin(X, axis=axis))

    else:
        mask = xp.isnan(X)
        X = xp.min(xp.where(mask, xp.asarray(+xp.inf, device=device(X)), X), axis=axis)
        # Replace Infs from all NaN slices with NaN again
        mask = xp.all(mask, axis=axis)
        if xp.any(mask):
            X = xp.where(mask, xp.asarray(xp.nan), X)
        return X


def _nanmax(X, axis=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _ = get_namespace(X)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmax(X, axis=axis))

    else:
        mask = xp.isnan(X)
        X = xp.max(xp.where(mask, xp.asarray(-xp.inf, device=device(X)), X), axis=axis)
        # Replace Infs from all NaN slices with NaN again
        mask = xp.all(mask, axis=axis)
        if xp.any(mask):
            X = xp.where(mask, xp.asarray(xp.nan), X)
        return X


def _asarray_with_order(array, dtype=None, order=None, copy=None, *, xp=None):
    """Helper to support the order kwarg only for NumPy-backed arrays

    Memory layout parameter `order` is not exposed in the Array API standard,
    however some input validation code in scikit-learn needs to work both
    for classes and functions that will leverage Array API only operations
    and for code that inherently relies on NumPy backed data containers with
    specific memory layout constraints (e.g. our own Cython code). The
    purpose of this helper is to make it possible to share code for data
    container validation without memory copies for both downstream use cases:
    the `order` parameter is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.
    """
    if xp is None:
        xp, _ = get_namespace(array)
    if _is_numpy_namespace(xp):
        # Use NumPy API to support order
        if copy is True:
            array = numpy.array(array, order=order, dtype=dtype)
        else:
            array = numpy.asarray(array, order=order, dtype=dtype)

        # At this point array is a NumPy ndarray. We convert it to an array
        # container that is consistent with the input's namespace.
        return xp.asarray(array)
    else:
        return xp.asarray(array, dtype=dtype, copy=copy)


def _convert_to_numpy(array, xp):
    """Convert X into a NumPy ndarray on the CPU."""
    xp_name = xp.__name__

    if xp_name in {"array_api_compat.torch", "torch"}:
        return array.cpu().numpy()
    elif xp_name == "cupy.array_api":
        return array._array.get()
    elif xp_name in {"array_api_compat.cupy", "cupy"}:  # pragma: nocover
        return array.get()

    return numpy.asarray(array)


def _estimator_with_converted_arrays(estimator, converter):
    """Create new estimator which converting all attributes that are arrays.

    The converter is called on all NumPy arrays and arrays that support the
    `DLPack interface <https://dmlc.github.io/dlpack/latest/>`__.

    Parameters
    ----------
    estimator : Estimator
        Estimator to convert

    converter : callable
        Callable that takes an array attribute and returns the converted array.

    Returns
    -------
    new_estimator : Estimator
        Convert estimator
    """
    from sklearn.base import clone

    new_estimator = clone(estimator)
    for key, attribute in vars(estimator).items():
        if hasattr(attribute, "__dlpack__") or isinstance(attribute, numpy.ndarray):
            attribute = converter(attribute)
        setattr(new_estimator, key, attribute)
    return new_estimator


def _atol_for_type(dtype):
    """Return the absolute tolerance for a given dtype."""
    return numpy.finfo(dtype).eps * 100
