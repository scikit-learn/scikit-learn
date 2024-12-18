"""Tools to support array_api."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import itertools
import math
import os
from functools import wraps

import numpy
import scipy
import scipy.sparse as sp
import scipy.special as special

from .._config import get_config
from .fixes import parse_version

_NUMPY_NAMESPACE_NAMES = {"numpy", "array_api_compat.numpy"}


def yield_namespaces(include_numpy_namespaces=True):
    """Yield supported namespace.

    This is meant to be used for testing purposes only.

    Parameters
    ----------
    include_numpy_namespaces : bool, default=True
        If True, also yield numpy namespaces.

    Returns
    -------
    array_namespace : str
        The name of the Array API namespace.
    """
    for array_namespace in [
        # The following is used to test the array_api_compat wrapper when
        # array_api_dispatch is enabled: in particular, the arrays used in the
        # tests are regular numpy arrays without any "device" attribute.
        "numpy",
        # Stricter NumPy-based Array API implementation. The
        # array_api_strict.Array instances always have a dummy "device" attribute.
        "array_api_strict",
        "cupy",
        "torch",
    ]:
        if not include_numpy_namespaces and array_namespace in _NUMPY_NAMESPACE_NAMES:
            continue
        yield array_namespace


def yield_namespace_device_dtype_combinations(include_numpy_namespaces=True):
    """Yield supported namespace, device, dtype tuples for testing.

    Use this to test that an estimator works with all combinations.

    Parameters
    ----------
    include_numpy_namespaces : bool, default=True
        If True, also yield numpy namespaces.

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
    for array_namespace in yield_namespaces(
        include_numpy_namespaces=include_numpy_namespaces
    ):
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
                f"NumPy must be {min_numpy_version} or newer (found"
                f" {numpy.__version__}) to dispatch array using"
                " the array API specification"
            )

        scipy_version = parse_version(scipy.__version__)
        min_scipy_version = "1.14.0"
        if scipy_version < parse_version(min_scipy_version):
            raise ImportError(
                f"SciPy must be {min_scipy_version} or newer"
                " (found {scipy.__version__}) to dispatch array using"
                " the array API specification"
            )

        if os.environ.get("SCIPY_ARRAY_API") != "1":
            raise RuntimeError(
                "Scikit-learn array API support was enabled but scipy's own support is "
                "not enabled. Please set the SCIPY_ARRAY_API=1 environment variable "
                "before importing sklearn or scipy. More details at: "
                "https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html"
            )


def _single_array_device(array):
    """Hardware device where the array data resides on."""
    if (
        isinstance(array, (numpy.ndarray, numpy.generic))
        or not hasattr(array, "device")
        # When array API dispatch is disabled, we expect the scikit-learn code
        # to use np.asarray so that the resulting NumPy array will implicitly use the
        # CPU. In this case, scikit-learn should stay as device neutral as possible,
        # hence the use of `device=None` which is accepted by all libraries, before
        # and after the expected conversion to NumPy via np.asarray.
        or not get_config()["array_api_dispatch"]
    ):
        return None
    else:
        return array.device


def device(*array_list, remove_none=True, remove_types=(str,)):
    """Hardware device where the array data resides on.

    If the hardware device is not the same for all arrays, an error is raised.

    Parameters
    ----------
    *array_list : arrays
        List of array instances from NumPy or an array API compatible library.

    remove_none : bool, default=True
        Whether to ignore None objects passed in array_list.

    remove_types : tuple or list, default=(str,)
        Types to ignore in array_list.

    Returns
    -------
    out : device
        `device` object (see the "Device Support" section of the array API spec).
    """
    array_list = _remove_non_arrays(
        *array_list, remove_none=remove_none, remove_types=remove_types
    )

    if not array_list:
        return None

    device_ = _single_array_device(array_list[0])

    # Note: here we cannot simply use a Python `set` as it requires
    # hashable members which is not guaranteed for Array API device
    # objects. In particular, CuPy devices are not hashable at the
    # time of writing.
    for array in array_list[1:]:
        device_other = _single_array_device(array)
        if device_ != device_other:
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
        # avoid circular import
        from ._unique import cached_unique

        a_unique, b_unique = cached_unique(a, b, xp=xp)
        return xp.asarray(numpy.union1d(a_unique, b_unique))
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
            return dtype in {xp.int8, xp.int16, xp.int32, xp.int64}
        elif kind == "unsigned integer":
            return dtype in {xp.uint8, xp.uint16, xp.uint32, xp.uint64}
        elif kind == "integral":
            return any(
                _isdtype_single(dtype, k, xp=xp)
                for k in ("signed integer", "unsigned integer")
            )
        elif kind == "real floating":
            return dtype in supported_float_dtypes(xp)
        elif kind == "complex floating":
            # Some name spaces might not have support for complex dtypes.
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


def supported_float_dtypes(xp):
    """Supported floating point types for the namespace.

    Note: float16 is not officially part of the Array API spec at the
    time of writing but scikit-learn estimators and functions can choose
    to accept it when xp.float16 is defined.

    https://data-apis.org/array-api/latest/API_specification/data_types.html
    """
    if hasattr(xp, "float16"):
        return (xp.float64, xp.float32, xp.float16)
    else:
        return (xp.float64, xp.float32)


def ensure_common_namespace_device(reference, *arrays):
    """Ensure that all arrays use the same namespace and device as reference.

    If necessary the arrays are moved to the same namespace and device as
    the reference array.

    Parameters
    ----------
    reference : array
        Reference array.

    *arrays : array
        Arrays to check.

    Returns
    -------
    arrays : list
        Arrays with the same namespace and device as reference.
    """
    xp, is_array_api = get_namespace(reference)

    if is_array_api:
        device_ = device(reference)
        # Move arrays to the same namespace and device as the reference array.
        return [xp.asarray(a, device=device_) for a in arrays]
    else:
        return arrays


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

    NumPy < 2 does not implement the namespace. NumPy 2 and later should
    progressively implement more an more of the latest Array API spec but this
    is still work in progress at this time.

    This wrapper makes it possible to write code that uses the standard Array
    API while working with any version of NumPy supported by scikit-learn.

    See the `get_namespace()` public function for more details.
    """

    # TODO: once scikit-learn drops support for NumPy < 2, this class can be
    # removed, assuming Array API compliance of NumPy 2 is actually sufficient
    # for scikit-learn's needs.

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

    def unique_all(self, x):
        return numpy.unique(
            x, return_index=True, return_inverse=True, return_counts=True
        )

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
        try:
            return isdtype(dtype, kind, xp=self)
        except TypeError:
            # In older versions of numpy, data types that arise from outside
            # numpy like from a Polars Series raise a TypeError.
            # e.g. TypeError: Cannot interpret 'Int64' as a data type.
            # Therefore, we return False.
            # TODO: Remove when minimum supported version of numpy is >= 1.21.
            return False

    def pow(self, x1, x2):
        return numpy.power(x1, x2)


_NUMPY_API_WRAPPER_INSTANCE = _NumPyAPIWrapper()


def _remove_non_arrays(*arrays, remove_none=True, remove_types=(str,)):
    """Filter arrays to exclude None and/or specific types.

    Raise ValueError if no arrays are left after filtering.

    Sparse arrays are always filtered out.

    Parameters
    ----------
    *arrays : array objects
        Array objects.

    remove_none : bool, default=True
        Whether to ignore None objects passed in arrays.

    remove_types : tuple or list, default=(str,)
        Types to ignore in the arrays.

    Returns
    -------
    filtered_arrays : list
        List of arrays filtered as requested. An empty list is returned if no input
        passes the filters.
    """
    filtered_arrays = []
    remove_types = tuple(remove_types)
    for array in arrays:
        if remove_none and array is None:
            continue
        if isinstance(array, remove_types):
            continue
        if sp.issparse(array):
            continue
        filtered_arrays.append(array)

    return filtered_arrays


def get_namespace(*arrays, remove_none=True, remove_types=(str,), xp=None):
    """Get namespace of arrays.

    Introspect `arrays` arguments and return their common Array API compatible
    namespace object, if any.

    Note that sparse arrays are filtered by default.

    See: https://numpy.org/neps/nep-0047-array-api-standard.html

    If `arrays` are regular numpy arrays, an instance of the `_NumPyAPIWrapper`
    compatibility wrapper is returned instead.

    Namespace support is not enabled by default. To enabled it call:

      sklearn.set_config(array_api_dispatch=True)

    or:

      with sklearn.config_context(array_api_dispatch=True):
          # your code here

    Otherwise an instance of the `_NumPyAPIWrapper` compatibility wrapper is
    always returned irrespective of the fact that arrays implement the
    `__array_namespace__` protocol or not.

    Note that if no arrays pass the set filters, ``_NUMPY_API_WRAPPER_INSTANCE, False``
    is returned.

    Parameters
    ----------
    *arrays : array objects
        Array objects.

    remove_none : bool, default=True
        Whether to ignore None objects passed in arrays.

    remove_types : tuple or list, default=(str,)
        Types to ignore in the arrays.

    xp : module, default=None
        Precomputed array namespace module. When passed, typically from a caller
        that has already performed inspection of its own inputs, skips array
        namespace inspection.

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
        if xp is not None:
            return xp, False
        else:
            return _NUMPY_API_WRAPPER_INSTANCE, False

    if xp is not None:
        return xp, True

    arrays = _remove_non_arrays(
        *arrays,
        remove_none=remove_none,
        remove_types=remove_types,
    )

    if not arrays:
        return _NUMPY_API_WRAPPER_INSTANCE, False

    _check_array_api_dispatch(array_api_dispatch)

    # array-api-compat is a required dependency of scikit-learn only when
    # configuring `array_api_dispatch=True`. Its import should therefore be
    # protected by _check_array_api_dispatch to display an informative error
    # message in case it is missing.
    import array_api_compat

    namespace, is_array_api_compliant = array_api_compat.get_namespace(*arrays), True

    if namespace.__name__ == "array_api_strict" and hasattr(
        namespace, "set_array_api_strict_flags"
    ):
        namespace.set_array_api_strict_flags(api_version="2023.12")

    return namespace, is_array_api_compliant


def get_namespace_and_device(*array_list, remove_none=True, remove_types=(str,)):
    """Combination into one single function of `get_namespace` and `device`.

    Parameters
    ----------
    *array_list : array objects
        Array objects.
    remove_none : bool, default=True
        Whether to ignore None objects passed in arrays.
    remove_types : tuple or list, default=(str,)
        Types to ignore in the arrays.

    Returns
    -------
    namespace : module
        Namespace shared by array objects. If any of the `arrays` are not arrays,
        the namespace defaults to NumPy.
    is_array_api_compliant : bool
        True if the arrays are containers that implement the Array API spec.
        Always False when array_api_dispatch=False.
    device : device
        `device` object (see the "Device Support" section of the array API spec).
    """
    array_list = _remove_non_arrays(
        *array_list,
        remove_none=remove_none,
        remove_types=remove_types,
    )

    skip_remove_kwargs = dict(remove_none=False, remove_types=[])

    xp, is_array_api = get_namespace(*array_list, **skip_remove_kwargs)
    arrays_device = device(*array_list, **skip_remove_kwargs)
    if is_array_api:
        return xp, is_array_api, arrays_device
    else:
        return xp, False, arrays_device


def _expit(X, xp=None):
    xp, _ = get_namespace(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(special.expit(numpy.asarray(X)))

    return 1.0 / (1.0 + xp.exp(-X))


def _fill_or_add_to_diagonal(array, value, xp, add_value=True, wrap=False):
    """Implementation to facilitate adding or assigning specified values to the
    diagonal of a 2-d array.

    If ``add_value`` is `True` then the values will be added to the diagonal
    elements otherwise the values will be assigned to the diagonal elements.
    By default, ``add_value`` is set to `True. This is currently only
    supported for 2-d arrays.

    The implementation is taken from the `numpy.fill_diagonal` function:
    https://github.com/numpy/numpy/blob/v2.0.0/numpy/lib/_index_tricks_impl.py#L799-L929
    """
    if array.ndim != 2:
        raise ValueError(
            f"array should be 2-d. Got array with shape {tuple(array.shape)}"
        )

    value = xp.asarray(value, dtype=array.dtype, device=device(array))
    end = None
    # Explicit, fast formula for the common case.  For 2-d arrays, we
    # accept rectangular ones.
    step = array.shape[1] + 1
    if not wrap:
        end = array.shape[1] * array.shape[1]

    array_flat = xp.reshape(array, (-1,))
    if add_value:
        array_flat[:end:step] += value
    else:
        array_flat[:end:step] = value


def _max_precision_float_dtype(xp, device):
    """Return the float dtype with the highest precision supported by the device."""
    # TODO: Update to use `__array_namespace__info__()` from array-api v2023.12
    # when/if that becomes more widespread.
    xp_name = xp.__name__
    if xp_name in {"array_api_compat.torch", "torch"} and (
        str(device).startswith("mps")
    ):  # pragma: no cover
        return xp.float32
    return xp.float64


def _find_matching_floating_dtype(*arrays, xp):
    """Find a suitable floating point dtype when computing with arrays.

    If any of the arrays are floating point, return the dtype with the highest
    precision by following official type promotion rules:

    https://data-apis.org/array-api/latest/API_specification/type_promotion.html

    If there are no floating point input arrays (all integral inputs for
    instance), return the default floating point dtype for the namespace.
    """
    dtyped_arrays = [a for a in arrays if hasattr(a, "dtype")]
    floating_dtypes = [
        a.dtype for a in dtyped_arrays if xp.isdtype(a.dtype, "real floating")
    ]
    if floating_dtypes:
        # Return the floating dtype with the highest precision:
        return xp.result_type(*floating_dtypes)

    # If none of the input arrays have a floating point dtype, they must be all
    # integer arrays or containers of Python scalars: return the default
    # floating point dtype for the namespace (implementation specific).
    return xp.asarray(0.0).dtype


def _average(a, axis=None, weights=None, normalize=True, xp=None):
    """Partial port of np.average to support the Array API.

    It does a best effort at mimicking the return dtype rule described at
    https://numpy.org/doc/stable/reference/generated/numpy.average.html but
    only for the common cases needed in scikit-learn.
    """
    xp, _, device_ = get_namespace_and_device(a, weights)

    if _is_numpy_namespace(xp):
        if normalize:
            return xp.asarray(numpy.average(a, axis=axis, weights=weights))
        elif axis is None and weights is not None:
            return xp.asarray(numpy.dot(a, weights))

    a = xp.asarray(a, device=device_)
    if weights is not None:
        weights = xp.asarray(weights, device=device_)

    if weights is not None and a.shape != weights.shape:
        if axis is None:
            raise TypeError(
                f"Axis must be specified when the shape of a {tuple(a.shape)} and "
                f"weights {tuple(weights.shape)} differ."
            )

        if tuple(weights.shape) != (a.shape[axis],):
            raise ValueError(
                f"Shape of weights weights.shape={tuple(weights.shape)} must be "
                f"consistent with a.shape={tuple(a.shape)} and {axis=}."
            )

        # If weights are 1D, add singleton dimensions for broadcasting
        shape = [1] * a.ndim
        shape[axis] = a.shape[axis]
        weights = xp.reshape(weights, shape)

    if xp.isdtype(a.dtype, "complex floating"):
        raise NotImplementedError(
            "Complex floating point values are not supported by average."
        )
    if weights is not None and xp.isdtype(weights.dtype, "complex floating"):
        raise NotImplementedError(
            "Complex floating point values are not supported by average."
        )

    output_dtype = _find_matching_floating_dtype(a, weights, xp=xp)
    a = xp.astype(a, output_dtype)

    if weights is None:
        return (xp.mean if normalize else xp.sum)(a, axis=axis)

    weights = xp.astype(weights, output_dtype)

    sum_ = xp.sum(xp.multiply(a, weights), axis=axis)

    if not normalize:
        return sum_

    scale = xp.sum(weights, axis=axis)
    if xp.any(scale == 0.0):
        raise ZeroDivisionError("Weights sum to zero, can't be normalized")

    return sum_ / scale


def _nanmin(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _ = get_namespace(X, xp=xp)
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


def _nanmax(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _ = get_namespace(X, xp=xp)
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


def _nanmean(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _ = get_namespace(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmean(X, axis=axis))
    else:
        mask = xp.isnan(X)
        total = xp.sum(xp.where(mask, xp.asarray(0.0, device=device(X)), X), axis=axis)
        count = xp.sum(xp.astype(xp.logical_not(mask), X.dtype), axis=axis)
        return total / count


def _asarray_with_order(
    array, dtype=None, order=None, copy=None, *, xp=None, device=None
):
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
    xp, _ = get_namespace(array, xp=xp)
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
        return xp.asarray(array, dtype=dtype, copy=copy, device=device)


def _ravel(array, xp=None):
    """Array API compliant version of np.ravel.

    For non numpy namespaces, it just returns a flattened array, that might
    be or not be a copy.
    """
    xp, _ = get_namespace(array, xp=xp)
    if _is_numpy_namespace(xp):
        array = numpy.asarray(array)
        return xp.asarray(numpy.ravel(array, order="C"))

    return xp.reshape(array, shape=(-1,))


def _convert_to_numpy(array, xp):
    """Convert X into a NumPy ndarray on the CPU."""
    xp_name = xp.__name__

    if xp_name in {"array_api_compat.torch", "torch"}:
        return array.cpu().numpy()
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


def _atol_for_type(dtype_or_dtype_name):
    """Return the absolute tolerance for a given numpy dtype."""
    if dtype_or_dtype_name is None:
        # If no dtype is specified when running tests for a given namespace, we
        # expect the same floating precision level as NumPy's default floating
        # point dtype.
        dtype_or_dtype_name = numpy.float64
    return numpy.finfo(dtype_or_dtype_name).eps * 100


def indexing_dtype(xp):
    """Return a platform-specific integer dtype suitable for indexing.

    On 32-bit platforms, this will typically return int32 and int64 otherwise.

    Note: using dtype is recommended for indexing transient array
    datastructures. For long-lived arrays, such as the fitted attributes of
    estimators, it is instead recommended to use platform-independent int32 if
    we do not expect to index more 2B elements. Using fixed dtypes simplifies
    the handling of serialized models, e.g. to deploy a model fit on a 64-bit
    platform to a target 32-bit platform such as WASM/pyodide.
    """
    # Currently this is implemented with simple hack that assumes that
    # following "may be" statements in the Array API spec always hold:
    # > The default integer data type should be the same across platforms, but
    # > the default may vary depending on whether Python is 32-bit or 64-bit.
    # > The default array index data type may be int32 on 32-bit platforms, but
    # > the default should be int64 otherwise.
    # https://data-apis.org/array-api/latest/API_specification/data_types.html#default-data-types
    # TODO: once sufficiently adopted, we might want to instead rely on the
    # newer inspection API: https://github.com/data-apis/array-api/issues/640
    return xp.asarray(0).dtype


def _searchsorted(a, v, *, side="left", sorter=None, xp=None):
    # Temporary workaround needed as long as searchsorted is not widely
    # adopted by implementers of the Array API spec. This is a quite
    # recent addition to the spec:
    # https://data-apis.org/array-api/latest/API_specification/generated/array_api.searchsorted.html # noqa
    xp, _ = get_namespace(a, v, xp=xp)
    if hasattr(xp, "searchsorted"):
        return xp.searchsorted(a, v, side=side, sorter=sorter)

    a_np = _convert_to_numpy(a, xp=xp)
    v_np = _convert_to_numpy(v, xp=xp)
    indices = numpy.searchsorted(a_np, v_np, side=side, sorter=sorter)
    return xp.asarray(indices, device=device(a))


def _setdiff1d(ar1, ar2, xp, assume_unique=False):
    """Find the set difference of two arrays.

    Return the unique values in `ar1` that are not in `ar2`.
    """
    if _is_numpy_namespace(xp):
        return xp.asarray(
            numpy.setdiff1d(
                ar1=ar1,
                ar2=ar2,
                assume_unique=assume_unique,
            )
        )

    if assume_unique:
        ar1 = xp.reshape(ar1, (-1,))
    else:
        ar1 = xp.unique_values(ar1)
        ar2 = xp.unique_values(ar2)
    return ar1[_in1d(ar1=ar1, ar2=ar2, xp=xp, assume_unique=True, invert=True)]


def _isin(element, test_elements, xp, assume_unique=False, invert=False):
    """Calculates ``element in test_elements``, broadcasting over `element`
    only.

    Returns a boolean array of the same shape as `element` that is True
    where an element of `element` is in `test_elements` and False otherwise.
    """
    if _is_numpy_namespace(xp):
        return xp.asarray(
            numpy.isin(
                element=element,
                test_elements=test_elements,
                assume_unique=assume_unique,
                invert=invert,
            )
        )

    original_element_shape = element.shape
    element = xp.reshape(element, (-1,))
    test_elements = xp.reshape(test_elements, (-1,))
    return xp.reshape(
        _in1d(
            ar1=element,
            ar2=test_elements,
            xp=xp,
            assume_unique=assume_unique,
            invert=invert,
        ),
        original_element_shape,
    )


# Note: This is a helper for the functions `_isin` and
# `_setdiff1d`. It is not meant to be called directly.
def _in1d(ar1, ar2, xp, assume_unique=False, invert=False):
    """Checks whether each element of an array is also present in a
    second array.

    Returns a boolean array the same length as `ar1` that is True
    where an element of `ar1` is in `ar2` and False otherwise.

    This function has been adapted using the original implementation
    present in numpy:
    https://github.com/numpy/numpy/blob/v1.26.0/numpy/lib/arraysetops.py#L524-L758
    """
    xp, _ = get_namespace(ar1, ar2, xp=xp)

    # This code is run to make the code significantly faster
    if ar2.shape[0] < 10 * ar1.shape[0] ** 0.145:
        if invert:
            mask = xp.ones(ar1.shape[0], dtype=xp.bool, device=device(ar1))
            for a in ar2:
                mask &= ar1 != a
        else:
            mask = xp.zeros(ar1.shape[0], dtype=xp.bool, device=device(ar1))
            for a in ar2:
                mask |= ar1 == a
        return mask

    if not assume_unique:
        ar1, rev_idx = xp.unique_inverse(ar1)
        ar2 = xp.unique_values(ar2)

    ar = xp.concat((ar1, ar2))
    device_ = device(ar)
    # We need this to be a stable sort.
    order = xp.argsort(ar, stable=True)
    reverse_order = xp.argsort(order, stable=True)
    sar = xp.take(ar, order, axis=0)
    if invert:
        bool_ar = sar[1:] != sar[:-1]
    else:
        bool_ar = sar[1:] == sar[:-1]
    flag = xp.concat((bool_ar, xp.asarray([invert], device=device_)))
    ret = xp.take(flag, reverse_order, axis=0)

    if assume_unique:
        return ret[: ar1.shape[0]]
    else:
        return xp.take(ret, rev_idx, axis=0)


def _count_nonzero(X, axis=None, sample_weight=None, xp=None, device=None):
    """A variant of `sklearn.utils.sparsefuncs.count_nonzero` for the Array API.

    If the array `X` is sparse, and we are using the numpy namespace then we
    simply call the original function. This function only supports 2D arrays.
    """
    from .sparsefuncs import count_nonzero

    xp, _ = get_namespace(X, sample_weight, xp=xp)
    if _is_numpy_namespace(xp) and sp.issparse(X):
        return count_nonzero(X, axis=axis, sample_weight=sample_weight)

    assert X.ndim == 2

    weights = xp.ones_like(X, device=device)
    if sample_weight is not None:
        sample_weight = xp.asarray(sample_weight, device=device)
        sample_weight = xp.reshape(sample_weight, (sample_weight.shape[0], 1))
        weights = xp.astype(weights, sample_weight.dtype) * sample_weight

    zero_scalar = xp.asarray(0, device=device, dtype=weights.dtype)
    return xp.sum(xp.where(X != 0, weights, zero_scalar), axis=axis)


def _modify_in_place_if_numpy(xp, func, *args, out=None, **kwargs):
    if _is_numpy_namespace(xp):
        func(*args, out=out, **kwargs)
    else:
        out = func(*args, **kwargs)
    return out


def _bincount(array, weights=None, minlength=None, xp=None):
    # TODO: update if bincount is ever adopted in a future version of the standard:
    # https://github.com/data-apis/array-api/issues/812
    xp, _ = get_namespace(array, xp=xp)
    if hasattr(xp, "bincount"):
        return xp.bincount(array, weights=weights, minlength=minlength)

    array_np = _convert_to_numpy(array, xp=xp)
    if weights is not None:
        weights_np = _convert_to_numpy(weights, xp=xp)
    else:
        weights_np = None
    bin_out = numpy.bincount(array_np, weights=weights_np, minlength=minlength)
    return xp.asarray(bin_out, device=device(array))


def _tolist(array, xp=None):
    xp, _ = get_namespace(array, xp=xp)
    if _is_numpy_namespace(xp):
        return array.tolist()
    array_np = _convert_to_numpy(array, xp=xp)
    return [element.item() for element in array_np]
