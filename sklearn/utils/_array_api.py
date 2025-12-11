"""Tools to support array_api."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import itertools
import math
import os

import numpy
import scipy
import scipy.sparse as sp
import scipy.special as special

from sklearn._config import get_config
from sklearn.externals import array_api_compat
from sklearn.externals import array_api_extra as xpx
from sklearn.externals.array_api_compat import numpy as np_compat
from sklearn.utils._dataframe import is_df_or_series
from sklearn.utils.fixes import parse_version

# TODO: complete __all__
__all__ = ["xpx"]  # we import xpx here just to re-export it, need this to appease ruff

_NUMPY_NAMESPACE_NAMES = {"numpy", "sklearn.externals.array_api_compat.numpy"}


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
    Use in conjunction with `ids=_get_namespace_device_dtype_ids` to give
    clearer pytest parametrization ID names.

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
                ("cpu", "cuda", "xpu"), ("float64", "float32")
            ):
                yield array_namespace, device, dtype
            yield array_namespace, "mps", "float32"

        elif array_namespace == "array_api_strict":
            try:
                import array_api_strict

                yield array_namespace, array_api_strict.Device("CPU_DEVICE"), "float64"
                yield array_namespace, array_api_strict.Device("device1"), "float32"
            except ImportError:
                # Those combinations will typically be skipped by pytest if
                # array_api_strict is not installed but we still need to see them in
                # the test output.
                yield array_namespace, "CPU_DEVICE", "float64"
                yield array_namespace, "device1", "float32"
        else:
            yield array_namespace, None, None


def _get_namespace_device_dtype_ids(param):
    """Get pytest parametrization IDs for `yield_namespace_device_dtype_combinations`"""
    # Gives clearer IDs for array-api-strict devices, see #31042 for details
    try:
        import array_api_strict
    except ImportError:
        # `None` results in the default pytest representation
        return None
    else:
        if param == array_api_strict.Device("CPU_DEVICE"):
            return "CPU_DEVICE"
        if param == array_api_strict.Device("device1"):
            return "device1"
        if param == array_api_strict.Device("device2"):
            return "device2"


def _check_array_api_dispatch(array_api_dispatch):
    """Checks that array API support is functional.

    In particular scipy needs to be recent enough and the environment variable
    needs to be set: SCIPY_ARRAY_API=1.
    """
    if not array_api_dispatch:
        return

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
        not hasattr(array, "device")
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
                f"Input arrays use different devices: {device_}, {device_other}"
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
        from sklearn.utils._unique import cached_unique

        a_unique, b_unique = cached_unique(a, b, xp=xp)
        return xp.asarray(numpy.union1d(a_unique, b_unique))
    assert a.ndim == b.ndim == 1
    return xp.unique_values(xp.concat([xp.unique_values(a), xp.unique_values(b)]))


def supported_float_dtypes(xp, device=None):
    """Supported floating point types for the namespace.

    Parameters
    ----------
    xp : module
        Array namespace to inspect.

    device : str or device instance from xp, default=None
        Device to use for dtype selection. If ``None``, then a default device
        is assumed.

    Returns
    -------
    supported_dtypes : tuple
        Tuple of real floating data types supported by the provided array namespace,
        ordered from the highest precision to lowest.

    See Also
    --------
    max_precision_float_dtype : Maximum float dtype for a namespace/device pair.

    Notes
    -----
    `float16` is not officially part of the Array API spec at the
    time of writing but scikit-learn estimators and functions can choose
    to accept it when xp.float16 is defined.

    Additionally, some devices available within a namespace may not support
    all floating-point types that the namespace provides.

    https://data-apis.org/array-api/latest/API_specification/data_types.html
    """
    dtypes_dict = xp.__array_namespace_info__().dtypes(
        kind="real floating", device=device
    )
    valid_float_dtypes = []
    for dtype_key in ("float64", "float32"):
        if dtype_key in dtypes_dict:
            valid_float_dtypes.append(dtypes_dict[dtype_key])

    if hasattr(xp, "float16"):
        valid_float_dtypes.append(xp.float16)

    return tuple(valid_float_dtypes)


def _remove_non_arrays(*arrays, remove_none=True, remove_types=(str,)):
    """Filter arrays to exclude None and/or specific types.

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
        if is_df_or_series(array):
            continue
        filtered_arrays.append(array)

    return filtered_arrays


def get_namespace(*arrays, remove_none=True, remove_types=(str,), xp=None):
    """Get namespace of arrays.

    Introspect `arrays` arguments and return their common Array API compatible
    namespace object, if any.

    Note that sparse arrays are filtered by default.

    See: https://numpy.org/neps/nep-0047-array-api-standard.html

    If `arrays` are regular numpy arrays, `array_api_compat.numpy` is returned instead.

    Namespace support is not enabled by default. To enabled it call:

      sklearn.set_config(array_api_dispatch=True)

    or:

      with sklearn.config_context(array_api_dispatch=True):
          # your code here

    Otherwise `array_api_compat.numpy` is
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
        the namespace defaults to the NumPy namespace.

    is_array_api_compliant : bool
        True if the arrays are containers that implement the array API spec (see
        https://data-apis.org/array-api/latest/index.html).
        Always False when array_api_dispatch=False.
    """
    array_api_dispatch = get_config()["array_api_dispatch"]
    if not array_api_dispatch:
        if xp is not None:
            return xp, False
        else:
            return np_compat, False

    if xp is not None:
        return xp, True

    arrays = _remove_non_arrays(
        *arrays,
        remove_none=remove_none,
        remove_types=remove_types,
    )

    if not arrays:
        return np_compat, False

    _check_array_api_dispatch(array_api_dispatch)

    namespace, is_array_api_compliant = array_api_compat.get_namespace(*arrays), True

    if namespace.__name__ == "array_api_strict" and hasattr(
        namespace, "set_array_api_strict_flags"
    ):
        namespace.set_array_api_strict_flags(api_version="2024.12")

    return namespace, is_array_api_compliant


def get_namespace_and_device(
    *array_list, remove_none=True, remove_types=(str,), xp=None
):
    """Combination into one single function of `get_namespace` and `device`.

    Parameters
    ----------
    *array_list : array objects
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
    device : device
        `device` object (see the "Device Support" section of the array API spec).
    """
    skip_remove_kwargs = dict(remove_none=False, remove_types=[])

    array_list = _remove_non_arrays(
        *array_list,
        remove_none=remove_none,
        remove_types=remove_types,
    )
    arrays_device = device(*array_list, **skip_remove_kwargs)

    if xp is None:
        xp, is_array_api = get_namespace(*array_list, **skip_remove_kwargs)
    else:
        xp, is_array_api = xp, True

    if is_array_api:
        return xp, is_array_api, arrays_device
    else:
        return xp, False, arrays_device


def move_to(*arrays, xp, device):
    """Move all arrays to `xp` and `device`.

    Each array will be moved to the reference namespace and device if
    it is not already using it. Otherwise the array is left unchanged.

    `array` may contain `None` entries, these are left unchanged.

    Sparse arrays are accepted (as pass through) if the reference namespace is
    NumPy, in which case they are returned unchanged. Otherwise a `TypeError`
    is raised.

    Parameters
    ----------
    *arrays : iterable of arrays
        Arrays to (potentially) move.

    xp : namespace
        Array API namespace to move arrays to.

    device : device
        Array API device to move arrays to.

    Returns
    -------
    arrays : tuple or array
        Tuple of arrays with the same namespace and device as reference. Single array
        returned if only one `arrays` input.
    """
    sparse_mask = [sp.issparse(array) for array in arrays]
    none_mask = [array is None for array in arrays]
    if any(sparse_mask) and not _is_numpy_namespace(xp):
        raise TypeError(
            "Sparse arrays are only accepted (and passed through) when the target "
            "namespace is Numpy"
        )

    converted_arrays = []

    for array, is_sparse, is_none in zip(arrays, sparse_mask, none_mask):
        if is_none:
            converted_arrays.append(None)
        elif is_sparse:
            converted_arrays.append(array)
        else:
            xp_array, _, device_array = get_namespace_and_device(array)
            if xp == xp_array and device == device_array:
                converted_arrays.append(array)
            else:
                try:
                    # The dlpack protocol is the future proof and library agnostic
                    # method to transfer arrays across namespace and device boundaries
                    # hence this method is attempted first and going through NumPy is
                    # only used as fallback in case of failure.
                    # Note: copy=None is the default since array-api 2023.12. Namespace
                    # libraries should only trigger a copy automatically if needed.
                    array_converted = xp.from_dlpack(array, device=device)
                    # `AttributeError` occurs when `__dlpack__` and `__dlpack_device__`
                    # methods are not present on the input array
                    # `TypeError` and `NotImplementedError` for packages that do not
                    # yet support dlpack 1.0
                    # (i.e. the `device`/`copy` kwargs, e.g., torch <= 2.8.0)
                    # See https://github.com/data-apis/array-api/pull/741 for
                    # more details about the introduction of the `copy` and `device`
                    # kwargs in the from_dlpack method and their expected
                    # meaning by namespaces implementing the array API spec.
                    # TODO: try removing this once DLPack v1 more widely supported
                    # TODO: ValueError should not be needed but is in practice:
                    # https://github.com/numpy/numpy/issues/30341
                except (
                    AttributeError,
                    TypeError,
                    NotImplementedError,
                    BufferError,
                    ValueError,
                ):
                    # Converting to numpy is tricky, handle this via dedicated function
                    if _is_numpy_namespace(xp):
                        array_converted = _convert_to_numpy(array, xp_array)
                    # Convert from numpy, all array libraries can do this
                    elif _is_numpy_namespace(xp_array):
                        array_converted = xp.asarray(array, device=device)
                    else:
                        # There is no generic way to convert from namespace A to B
                        # So we first convert from A to numpy and then from numpy to B
                        # The way to avoid this round trip is to lobby for DLpack
                        # support in libraries A and B
                        array_np = _convert_to_numpy(array, xp_array)
                        array_converted = xp.asarray(array_np, device=device)
                converted_arrays.append(array_converted)

    return (
        converted_arrays[0] if len(converted_arrays) == 1 else tuple(converted_arrays)
    )


def _expit(X, xp=None):
    xp, _ = get_namespace(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(special.expit(numpy.asarray(X)))

    return 1.0 / (1.0 + xp.exp(-X))


def _validate_diagonal_args(array, value, xp):
    """Validate arguments to `_fill_diagonal`/`_add_to_diagonal`."""
    if array.ndim != 2:
        raise ValueError(
            f"`array` should be 2D. Got array with shape {tuple(array.shape)}"
        )

    value = xp.asarray(value, dtype=array.dtype, device=device(array))
    if value.ndim not in [0, 1]:
        raise ValueError(
            "`value` needs to be a scalar or a 1D array, "
            f"got a {value.ndim}D array instead."
        )
    min_rows_columns = min(array.shape)
    if value.ndim == 1 and value.shape[0] != min_rows_columns:
        raise ValueError(
            "`value` needs to be a scalar or 1D array of the same length as the "
            f"diagonal of `array` ({min_rows_columns}). Got {value.shape[0]}"
        )

    return value, min_rows_columns


def _fill_diagonal(array, value, xp):
    """Minimal implementation of `numpy.fill_diagonal`.

    `wrap` is not supported (i.e. always False). `value` should be a scalar or
    1D of greater or equal length as the diagonal (i.e., `value` is never repeated
    when shorter).

    Note `array` is altered in place.
    """
    value, min_rows_columns = _validate_diagonal_args(array, value, xp)

    if _is_numpy_namespace(xp):
        xp.fill_diagonal(array, value, wrap=False)
    else:
        # TODO: when array libraries support `reshape(copy)`, use
        # `reshape(array, (-1,), copy=False)`, then fill with `[:end:step]` (within
        # `try/except`). This is faster than for loop, when no copy needs to be
        # made within `reshape`. See #31445 for details.
        if value.ndim == 0:
            for i in range(min_rows_columns):
                array[i, i] = value
        else:
            for i in range(min_rows_columns):
                array[i, i] = value[i]


def _add_to_diagonal(array, value, xp):
    """Add `value` to diagonal of `array`.

    Related to `fill_diagonal`. `value` should be a scalar or
    1D of greater or equal length as the diagonal (i.e., `value` is never repeated
    when shorter).

    Note `array` is altered in place.
    """
    value, min_rows_columns = _validate_diagonal_args(array, value, xp)

    if _is_numpy_namespace(xp):
        step = array.shape[1] + 1
        # Ensure we do not wrap
        end = array.shape[1] * array.shape[1]
        array.flat[:end:step] += value
        return

    # TODO: when array libraries support `reshape(copy)`, use
    # `reshape(array, (-1,), copy=False)`, then fill with `[:end:step]` (within
    # `try/except`). This is faster than for loop, when no copy needs to be
    # made within `reshape`. See #31445 for details.
    value = xp.linalg.diagonal(array) + value
    for i in range(min_rows_columns):
        array[i, i] = value[i]


def _is_xp_namespace(xp, name):
    return xp.__name__ in (
        name,
        f"array_api_compat.{name}",
        f"sklearn.externals.array_api_compat.{name}",
    )


def _max_precision_float_dtype(xp, device):
    """Return the float dtype with the highest precision supported by the device."""
    # TODO: Update to use `__array_namespace__info__()` from array-api v2023.12
    # when/if that becomes more widespread.
    if _is_xp_namespace(xp, "torch") and str(device).startswith(
        "mps"
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
    dtyped_arrays = [xp.asarray(a) for a in arrays if hasattr(a, "dtype")]
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
    xp, _, device_ = get_namespace_and_device(a, weights, xp=xp)

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
        weights = xp.reshape(weights, tuple(shape))

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


def _median(x, axis=None, keepdims=False, xp=None):
    # XXX: `median` is not included in the array API spec, but is implemented
    # in most array libraries, and all that we support (as of May 2025).
    # TODO: consider simplifying this code to use scipy instead once the oldest
    # supported SciPy version provides `scipy.stats.quantile` with native array API
    # support (likely scipy 1.16 at the time of writing). Proper benchmarking of
    # either option with popular array namespaces is required to evaluate the
    # impact of this choice.
    xp, _, device = get_namespace_and_device(x, xp=xp)

    # `torch.median` takes the lower of the two medians when `x` has even number
    # of elements, thus we use `torch.quantile(q=0.5)`, which gives mean of the two
    if array_api_compat.is_torch_namespace(xp):
        return xp.quantile(x, q=0.5, dim=axis, keepdim=keepdims)

    if hasattr(xp, "median"):
        return xp.median(x, axis=axis, keepdims=keepdims)

    # Intended mostly for array-api-strict (which as no "median", as per the spec)
    # as `_convert_to_numpy` does not necessarily work for all array types.
    x_np = _convert_to_numpy(x, xp=xp)
    return xp.asarray(numpy.median(x_np, axis=axis, keepdims=keepdims), device=device)


def _xlogy(x, y, xp=None):
    # TODO: Remove this once https://github.com/scipy/scipy/issues/21736 is fixed
    xp, _, device_ = get_namespace_and_device(x, y, xp=xp)

    with numpy.errstate(divide="ignore", invalid="ignore"):
        temp = x * xp.log(y)
    return xp.where(x == 0.0, xp.asarray(0.0, dtype=temp.dtype, device=device_), temp)


def _nanmin(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _, device_ = get_namespace_and_device(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmin(X, axis=axis))

    else:
        mask = xp.isnan(X)
        X = xp.min(
            xp.where(mask, xp.asarray(+xp.inf, dtype=X.dtype, device=device_), X),
            axis=axis,
        )
        # Replace Infs from all NaN slices with NaN again
        mask = xp.all(mask, axis=axis)
        if xp.any(mask):
            X = xp.where(mask, xp.asarray(xp.nan, dtype=X.dtype, device=device_), X)
        return X


def _nanmax(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _, device_ = get_namespace_and_device(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmax(X, axis=axis))

    else:
        mask = xp.isnan(X)
        X = xp.max(
            xp.where(mask, xp.asarray(-xp.inf, dtype=X.dtype, device=device_), X),
            axis=axis,
        )
        # Replace Infs from all NaN slices with NaN again
        mask = xp.all(mask, axis=axis)
        if xp.any(mask):
            X = xp.where(mask, xp.asarray(xp.nan, dtype=X.dtype, device=device_), X)
        return X


def _nanmean(X, axis=None, xp=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _, device_ = get_namespace_and_device(X, xp=xp)
    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nanmean(X, axis=axis))
    else:
        mask = xp.isnan(X)
        total = xp.sum(
            xp.where(mask, xp.asarray(0.0, dtype=X.dtype, device=device_), X), axis=axis
        )
        count = xp.sum(xp.astype(xp.logical_not(mask), X.dtype), axis=axis)
        return total / count


def _nansum(X, axis=None, xp=None, keepdims=False, dtype=None):
    # TODO: refactor once nan-aware reductions are standardized:
    # https://github.com/data-apis/array-api/issues/621
    xp, _, X_device = get_namespace_and_device(X, xp=xp)

    if _is_numpy_namespace(xp):
        return xp.asarray(numpy.nansum(X, axis=axis, keepdims=keepdims, dtype=dtype))

    mask = xp.isnan(X)
    masked_arr = xp.where(mask, xp.asarray(0, device=X_device, dtype=X.dtype), X)
    return xp.sum(masked_arr, axis=axis, keepdims=keepdims, dtype=dtype)


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
    if _is_xp_namespace(xp, "torch"):
        return array.cpu().numpy()
    elif _is_xp_namespace(xp, "cupy"):  # pragma: nocover
        return array.get()
    elif _is_xp_namespace(xp, "array_api_strict"):
        return numpy.asarray(xp.asarray(array, device=xp.Device("CPU_DEVICE")))

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
    return numpy.finfo(dtype_or_dtype_name).eps * 1000


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


# Note: This is a helper for the function `_isin`.
# It is not meant to be called directly.
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
    if size(sar) >= 1:
        bool_ar = sar[1:] != sar[:-1] if invert else sar[1:] == sar[:-1]
    else:
        # indexing undefined in standard when sar is empty
        bool_ar = xp.asarray([False]) if invert else xp.asarray([True])
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
    from sklearn.utils.sparsefuncs import count_nonzero

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


def _logsumexp(array, axis=None, xp=None):
    # TODO replace by scipy.special.logsumexp when
    # https://github.com/scipy/scipy/pull/22683 is part of a release.
    # The following code is strongly inspired and simplified from
    # scipy.special._logsumexp.logsumexp
    xp, _, device = get_namespace_and_device(array, xp=xp)
    axis = tuple(range(array.ndim)) if axis is None else axis

    supported_dtypes = supported_float_dtypes(xp)
    if array.dtype not in supported_dtypes:
        array = xp.asarray(array, dtype=supported_dtypes[0])

    array_max = xp.max(array, axis=axis, keepdims=True)
    index_max = array == array_max

    array = xp.asarray(array, copy=True)
    array[index_max] = -xp.inf
    i_max_dt = xp.astype(index_max, array.dtype)
    m = xp.sum(i_max_dt, axis=axis, keepdims=True, dtype=array.dtype)
    # Specifying device explicitly is the fix for https://github.com/scipy/scipy/issues/22680
    shift = xp.where(
        xp.isfinite(array_max),
        array_max,
        xp.asarray(0, dtype=array_max.dtype, device=device),
    )
    exp = xp.exp(array - shift)
    s = xp.sum(exp, axis=axis, keepdims=True, dtype=exp.dtype)
    s = xp.where(s == 0, s, s / m)
    out = xp.log1p(s) + xp.log(m) + array_max
    out = xp.squeeze(out, axis=axis)
    out = out[()] if out.ndim == 0 else out

    return out


def _cholesky(covariance, xp):
    if _is_numpy_namespace(xp):
        return scipy.linalg.cholesky(covariance, lower=True)
    else:
        return xp.linalg.cholesky(covariance)


def _linalg_solve(cov_chol, eye_matrix, xp):
    if _is_numpy_namespace(xp):
        return scipy.linalg.solve_triangular(cov_chol, eye_matrix, lower=True)
    else:
        return xp.linalg.solve(cov_chol, eye_matrix)


def _half_multinomial_loss(y, pred, sample_weight=None, xp=None):
    """A version of the multinomial loss that is compatible with the array API"""
    xp, _, device_ = get_namespace_and_device(y, pred, sample_weight)
    log_sum_exp = _logsumexp(pred, axis=1, xp=xp)
    y = xp.asarray(y, dtype=xp.int64, device=device_)
    class_margins = xp.arange(y.shape[0], device=device_) * pred.shape[1]
    label_predictions = xp.take(_ravel(pred), y + class_margins)
    return float(
        _average(log_sum_exp - label_predictions, weights=sample_weight, xp=xp)
    )
