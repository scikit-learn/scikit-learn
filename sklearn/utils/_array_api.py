"""Tools to support array_api."""
from functools import wraps
import math

import numpy
import scipy.special as special

from .._config import get_config


def device(x):
    """
    Hardware device the array data resides on.

    Parameters
    ----------
    x: array
        array instance from NumPy or an array API compatible library.

    Returns
    -------
    out: device
        a ``device`` object (see the "Device Support" section of the array API spec).
    """
    if isinstance(x, (numpy.ndarray, numpy.generic)):
        return "cpu"
    return x.device


def size(x):
    """
    Return the total number of elements of x
    """
    if None in x.shape:
        return None
    return math.prod(x.shape)


def _is_numpy_namespace(xp):
    return xp.__name__ in {
        "numpy",
        "_NumPyApiWrapper",
        "array_api_compat.numpy",
        "numpy.array_api",
    }


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
            return isdtype(dtype, ("signed integer", "unsigned integer"), xp=xp)
        elif kind == "real floating":
            return dtype in {xp.float32, xp.float64}
        elif kind == "complex floating":
            # Some name spaces do not have complex, such as cupy.array_api
            # and numpy.array_api
            if hasattr(xp, "complex64"):
                return dtype == xp.complex64
            if hasattr(xp, "complex128"):
                return dtype == xp.complex128
            return False
        elif kind == "numeric":
            return isdtype(
                dtype, ("integral", "real floating", "complex floating"), xp=xp
            )
        else:
            raise ValueError(f"Unrecognized data type kind: {kind!r}")
    else:
        return dtype == kind


class _ArrayAPIWrapper:
    """sklearn specific Array API compatibility wrapper

    This wrapper makes it possible for scikit-learn maintainers to
    deal with discrepancies between different implementations of the
    Python array API standard and its evolution over time.

    The Python array API standard specification:
    https://data-apis.org/array-api/latest/

    Documentation of the NumPy implementation:
    https://numpy.org/neps/nep-0047-array-api-standard.html
    """

    def __init__(self, array_namespace):
        self._namespace = array_namespace

    def __getattr__(self, name):
        return getattr(self._namespace, name)

    def take(self, X, indices, *, axis=0):
        # When array_api supports `take` we can use this directly
        # https://github.com/data-apis/array-api/issues/177
        if self._namespace.__name__ == "numpy.array_api":
            X_np = numpy.take(X, indices, axis=axis)
            return self._namespace.asarray(X_np)

        # We only support axis in (0, 1) and ndim in (1, 2) because that is all we need
        # in scikit-learn
        if axis not in {0, 1}:
            raise ValueError(f"Only axis in (0, 1) is supported. Got {axis}")

        if X.ndim not in {1, 2}:
            raise ValueError(f"Only X.ndim in (1, 2) is supported. Got {X.ndim}")

        if axis == 0:
            if X.ndim == 1:
                selected = [X[i] for i in indices]
            else:  # X.ndim == 2
                selected = [X[i, :] for i in indices]
        else:  # axis == 1
            selected = [X[:, i] for i in indices]
        return self._namespace.stack(selected, axis=axis)

    def isdtype(self, dtype, kind):
        """Returns a boolean indicating whether a provided dtype is of type "kind".

        Included in the v2022.12 of the Array API spec.
        https://data-apis.org/array-api/latest/API_specification/generated/array_api.isdtype.html
        """
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


class _NumPyApiWrapper:
    """Array API compat wrapper for any numpy version

    NumPy < 1.22 does not expose the numpy.array_api namespace. This
    wrapper makes it possible to write code that uses the standard
    Array API while working with any version of NumPy supported by
    scikit-learn.

    See the `get_namespace()` public function for more details.
    """

    def __getattr__(self, name):
        attr = getattr(numpy, name)
        # accept device
        if name in {
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
        }:
            return _accept_device_cpu(attr)

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

    def isdtype(self, dtype, kind):
        """Returns a boolean indicating whether a provided dtype is of type "kind".

        Included in the v2022.12 of the Array API spec.
        https://data-apis.org/array-api/latest/API_specification/generated/array_api.isdtype.html
        """
        return isdtype(dtype, kind, xp=self)


_NUMPY_API_WRAPPER_INSTANCE = _NumPyApiWrapper()


def get_namespace(*arrays):
    """Get namespace of arrays.

    Introspect `arrays` arguments and return their common Array API
    compatible namespace object, if any. NumPy 1.22 and later can
    construct such containers using the `numpy.array_api` namespace
    for instance.

    See: https://numpy.org/neps/nep-0047-array-api-standard.html

    If `arrays` are regular numpy arrays, an instance of the
    `_NumPyApiWrapper` compatibility wrapper is returned instead.

    Namespace support is not enabled by default. To enabled it
    call:

      sklearn.set_config(array_api_dispatch=True)

    or:

      with sklearn.config_context(array_api_dispatch=True):
          # your code here

    Otherwise an instance of the `_NumPyApiWrapper`
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

    is_array_api : bool
        True of the arrays are containers that implement the Array API spec.
    """
    array_api_dispatch = get_config()["array_api_dispatch"]
    return _get_namespace(*arrays, array_api_dispatch=array_api_dispatch)


def _get_namespace(*arrays, array_api_dispatch):
    """Helper method for get_namespace that dispatches with array_api_dispatch.

    Parameters
    ----------
    *arrays : array objects
        Array objects.

    array_api_dispatch : bool
        If True, the array namespace is obtained from the array objects.

    Returns
    -------
    namespace : module
        Namespace shared by array objects. If any of the `arrays` are not arrays,
        the namespace defaults to NumPy.

    is_array_api : bool
        True of the arrays are containers that implement the Array API spec.
    """
    if not array_api_dispatch:
        return _NUMPY_API_WRAPPER_INSTANCE, False

    try:
        import array_api_compat
    except ImportError:
        return _NUMPY_API_WRAPPER_INSTANCE, False

    try:
        namespace, is_array_api = array_api_compat.get_namespace(*arrays), True
    except TypeError:
        return _NUMPY_API_WRAPPER_INSTANCE, False

    if namespace.__name__ in {"numpy.array_api", "cupy.array_api"}:
        namespace = _ArrayAPIWrapper(namespace)

    return namespace, is_array_api


def _expit(X):
    xp, _ = get_namespace(X)
    if _is_numpy_namespace(xp):
        return xp.asarray(special.expit(numpy.asarray(X)))

    return 1.0 / (1.0 + xp.exp(-X))


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
        return xp.asarray(array)
    else:
        return xp.asarray(array, dtype=dtype, copy=copy)


def _convert_to_numpy(array):
    """Convert X into a NumPy ndarray on the CPU."""
    xp, _ = _get_namespace(array, array_api_dispatch=True)

    xp_name = xp.__name__

    if xp_name in {"array_api_compat.torch", "torch"}:
        return array.cpu().numpy()
    elif xp_name == "cupy.array_api":
        return array._array.get()
    elif xp_name in {"array_api_compat.cupy", "cupy"}:
        return array.get()

    return numpy.asarray(array)


def _estimator_with_converted_arrays(estimator, converter):
    """Create new estimator which converting all attributes that are arrays.

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
        _, is_array_api = _get_namespace(attribute, array_api_dispatch=True)
        if is_array_api:
            attribute = converter(attribute)
        setattr(new_estimator, key, attribute)
    return new_estimator
