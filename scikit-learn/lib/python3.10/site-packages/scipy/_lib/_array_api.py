"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
import os

from types import ModuleType
from typing import Any, Literal, TypeAlias

import numpy as np
import numpy.typing as npt

from scipy._lib import array_api_compat
from scipy._lib.array_api_compat import (
    is_array_api_obj,
    size as xp_size,
    numpy as np_compat,
    device as xp_device,
    is_numpy_namespace as is_numpy,
    is_cupy_namespace as is_cupy,
    is_torch_namespace as is_torch,
    is_jax_namespace as is_jax,
    is_array_api_strict_namespace as is_array_api_strict
)

__all__ = [
    '_asarray', 'array_namespace', 'assert_almost_equal', 'assert_array_almost_equal',
    'get_xp_devices',
    'is_array_api_strict', 'is_complex', 'is_cupy', 'is_jax', 'is_numpy', 'is_torch', 
    'SCIPY_ARRAY_API', 'SCIPY_DEVICE', 'scipy_namespace_for',
    'xp_assert_close', 'xp_assert_equal', 'xp_assert_less',
    'xp_copy', 'xp_copysign', 'xp_device',
    'xp_moveaxis_to_end', 'xp_ravel', 'xp_real', 'xp_sign', 'xp_size',
    'xp_take_along_axis', 'xp_unsupported_param_msg', 'xp_vector_norm',
]


# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")

_GLOBAL_CONFIG = {
    "SCIPY_ARRAY_API": SCIPY_ARRAY_API,
    "SCIPY_DEVICE": SCIPY_DEVICE,
}


Array: TypeAlias = Any  # To be changed to a Protocol later (see array-api#589)
ArrayLike: TypeAlias = Array | npt.ArrayLike


def _compliance_scipy(arrays):
    """Raise exceptions on known-bad subclasses.

    The following subclasses are not supported and raise and error:
    - `numpy.ma.MaskedArray`
    - `numpy.matrix`
    - NumPy arrays which do not have a boolean or numerical dtype
    - Any array-like which is neither array API compatible nor coercible by NumPy
    - Any array-like which is coerced by NumPy to an unsupported dtype
    """
    for i in range(len(arrays)):
        array = arrays[i]

        from scipy.sparse import issparse
        # this comes from `_util._asarray_validated`
        if issparse(array):
            msg = ('Sparse arrays/matrices are not supported by this function. '
                   'Perhaps one of the `scipy.sparse.linalg` functions '
                   'would work instead.')
            raise ValueError(msg)

        if isinstance(array, np.ma.MaskedArray):
            raise TypeError("Inputs of type `numpy.ma.MaskedArray` are not supported.")
        elif isinstance(array, np.matrix):
            raise TypeError("Inputs of type `numpy.matrix` are not supported.")
        if isinstance(array, np.ndarray | np.generic):
            dtype = array.dtype
            if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
                raise TypeError(f"An argument has dtype `{dtype!r}`; "
                                f"only boolean and numerical dtypes are supported.")
        elif not is_array_api_obj(array):
            try:
                array = np.asanyarray(array)
            except TypeError:
                raise TypeError("An argument is neither array API compatible nor "
                                "coercible by NumPy.")
            dtype = array.dtype
            if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
                message = (
                    f"An argument was coerced to an unsupported dtype `{dtype!r}`; "
                    f"only boolean and numerical dtypes are supported."
                )
                raise TypeError(message)
            arrays[i] = array
    return arrays


def _check_finite(array: Array, xp: ModuleType) -> None:
    """Check for NaNs or Infs."""
    msg = "array must not contain infs or NaNs"
    try:
        if not xp.all(xp.isfinite(array)):
            raise ValueError(msg)
    except TypeError:
        raise ValueError(msg)


def array_namespace(*arrays: Array) -> ModuleType:
    """Get the array API compatible namespace for the arrays xs.

    Parameters
    ----------
    *arrays : sequence of array_like
        Arrays used to infer the common namespace.

    Returns
    -------
    namespace : module
        Common namespace.

    Notes
    -----
    Thin wrapper around `array_api_compat.array_namespace`.

    1. Check for the global switch: SCIPY_ARRAY_API. This can also be accessed
       dynamically through ``_GLOBAL_CONFIG['SCIPY_ARRAY_API']``.
    2. `_compliance_scipy` raise exceptions on known-bad subclasses. See
       its definition for more details.

    When the global switch is False, it defaults to the `numpy` namespace.
    In that case, there is no compliance check. This is a convenience to
    ease the adoption. Otherwise, arrays must comply with the new rules.
    """
    if not _GLOBAL_CONFIG["SCIPY_ARRAY_API"]:
        # here we could wrap the namespace if needed
        return np_compat

    _arrays = [array for array in arrays if array is not None]

    _arrays = _compliance_scipy(_arrays)

    return array_api_compat.array_namespace(*_arrays)


def _asarray(
        array: ArrayLike,
        dtype: Any = None,
        order: Literal['K', 'A', 'C', 'F'] | None = None,
        copy: bool | None = None,
        *,
        xp: ModuleType | None = None,
        check_finite: bool = False,
        subok: bool = False,
    ) -> Array:
    """SciPy-specific replacement for `np.asarray` with `order`, `check_finite`, and
    `subok`.

    Memory layout parameter `order` is not exposed in the Array API standard.
    `order` is only enforced if the input array implementation
    is NumPy based, otherwise `order` is just silently ignored.

    `check_finite` is also not a keyword in the array API standard; included
    here for convenience rather than that having to be a separate function
    call inside SciPy functions.

    `subok` is included to allow this function to preserve the behaviour of
    `np.asanyarray` for NumPy based inputs.
    """
    if xp is None:
        xp = array_namespace(array)
    if is_numpy(xp):
        # Use NumPy API to support order
        if copy is True:
            array = np.array(array, order=order, dtype=dtype, subok=subok)
        elif subok:
            array = np.asanyarray(array, order=order, dtype=dtype)
        else:
            array = np.asarray(array, order=order, dtype=dtype)
    else:
        try:
            array = xp.asarray(array, dtype=dtype, copy=copy)
        except TypeError:
            coerced_xp = array_namespace(xp.asarray(3))
            array = coerced_xp.asarray(array, dtype=dtype, copy=copy)

    if check_finite:
        _check_finite(array, xp)

    return array


def xp_copy(x: Array, *, xp: ModuleType | None = None) -> Array:
    """
    Copies an array.

    Parameters
    ----------
    x : array

    xp : array_namespace

    Returns
    -------
    copy : array
        Copied array

    Notes
    -----
    This copy function does not offer all the semantics of `np.copy`, i.e. the
    `subok` and `order` keywords are not used.
    """
    # Note: for older NumPy versions, `np.asarray` did not support the `copy` kwarg,
    # so this uses our other helper `_asarray`.
    if xp is None:
        xp = array_namespace(x)

    return _asarray(x, copy=True, xp=xp)


def _strict_check(actual, desired, xp, *,
                  check_namespace=True, check_dtype=True, check_shape=True,
                  check_0d=True):
    __tracebackhide__ = True  # Hide traceback for py.test
    if check_namespace:
        _assert_matching_namespace(actual, desired)

    # only NumPy distinguishes between scalars and arrays; we do if check_0d=True.
    # do this first so we can then cast to array (and thus use the array API) below.
    if is_numpy(xp) and check_0d:
        _msg = ("Array-ness does not match:\n Actual: "
                f"{type(actual)}\n Desired: {type(desired)}")
        assert ((xp.isscalar(actual) and xp.isscalar(desired))
                or (not xp.isscalar(actual) and not xp.isscalar(desired))), _msg

    actual = xp.asarray(actual)
    desired = xp.asarray(desired)

    if check_dtype:
        _msg = f"dtypes do not match.\nActual: {actual.dtype}\nDesired: {desired.dtype}"
        assert actual.dtype == desired.dtype, _msg

    if check_shape:
        _msg = f"Shapes do not match.\nActual: {actual.shape}\nDesired: {desired.shape}"
        assert actual.shape == desired.shape, _msg

    desired = xp.broadcast_to(desired, actual.shape)
    return actual, desired


def _assert_matching_namespace(actual, desired):
    __tracebackhide__ = True  # Hide traceback for py.test
    actual = actual if isinstance(actual, tuple) else (actual,)
    desired_space = array_namespace(desired)
    for arr in actual:
        arr_space = array_namespace(arr)
        _msg = (f"Namespaces do not match.\n"
                f"Actual: {arr_space.__name__}\n"
                f"Desired: {desired_space.__name__}")
        assert arr_space == desired_space, _msg


def xp_assert_equal(actual, desired, *, check_namespace=True, check_dtype=True,
                    check_shape=True, check_0d=True, err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)

    actual, desired = _strict_check(
        actual, desired, xp, check_namespace=check_namespace,
        check_dtype=check_dtype, check_shape=check_shape,
        check_0d=check_0d
    )

    if is_cupy(xp):
        return xp.testing.assert_array_equal(actual, desired, err_msg=err_msg)
    elif is_torch(xp):
        # PyTorch recommends using `rtol=0, atol=0` like this
        # to test for exact equality
        err_msg = None if err_msg == '' else err_msg
        return xp.testing.assert_close(actual, desired, rtol=0, atol=0, equal_nan=True,
                                       check_dtype=False, msg=err_msg)
    # JAX uses `np.testing`
    return np.testing.assert_array_equal(actual, desired, err_msg=err_msg)


def xp_assert_close(actual, desired, *, rtol=None, atol=0, check_namespace=True,
                    check_dtype=True, check_shape=True, check_0d=True,
                    err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)

    actual, desired = _strict_check(
        actual, desired, xp,
        check_namespace=check_namespace, check_dtype=check_dtype,
        check_shape=check_shape, check_0d=check_0d
    )

    floating = xp.isdtype(actual.dtype, ('real floating', 'complex floating'))
    if rtol is None and floating:
        # multiplier of 4 is used as for `np.float64` this puts the default `rtol`
        # roughly half way between sqrt(eps) and the default for
        # `numpy.testing.assert_allclose`, 1e-7
        rtol = xp.finfo(actual.dtype).eps**0.5 * 4
    elif rtol is None:
        rtol = 1e-7

    if is_cupy(xp):
        return xp.testing.assert_allclose(actual, desired, rtol=rtol,
                                          atol=atol, err_msg=err_msg)
    elif is_torch(xp):
        err_msg = None if err_msg == '' else err_msg
        return xp.testing.assert_close(actual, desired, rtol=rtol, atol=atol,
                                       equal_nan=True, check_dtype=False, msg=err_msg)
    # JAX uses `np.testing`
    return np.testing.assert_allclose(actual, desired, rtol=rtol,
                                      atol=atol, err_msg=err_msg)


def xp_assert_less(actual, desired, *, check_namespace=True, check_dtype=True,
                   check_shape=True, check_0d=True, err_msg='', verbose=True, xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)

    actual, desired = _strict_check(
        actual, desired, xp, check_namespace=check_namespace,
        check_dtype=check_dtype, check_shape=check_shape,
        check_0d=check_0d
    )

    if is_cupy(xp):
        return xp.testing.assert_array_less(actual, desired,
                                            err_msg=err_msg, verbose=verbose)
    elif is_torch(xp):
        if actual.device.type != 'cpu':
            actual = actual.cpu()
        if desired.device.type != 'cpu':
            desired = desired.cpu()
    # JAX uses `np.testing`
    return np.testing.assert_array_less(actual, desired,
                                        err_msg=err_msg, verbose=verbose)


def assert_array_almost_equal(actual, desired, decimal=6, *args, **kwds):
    """Backwards compatible replacement. In new code, use xp_assert_close instead.
    """
    rtol, atol = 0, 1.5*10**(-decimal)
    return xp_assert_close(actual, desired,
                           atol=atol, rtol=rtol, check_dtype=False, check_shape=False,
                           *args, **kwds)


def assert_almost_equal(actual, desired, decimal=7, *args, **kwds):
    """Backwards compatible replacement. In new code, use xp_assert_close instead.
    """
    rtol, atol = 0, 1.5*10**(-decimal)
    return xp_assert_close(actual, desired,
                           atol=atol, rtol=rtol, check_dtype=False, check_shape=False,
                           *args, **kwds)


def xp_unsupported_param_msg(param: Any) -> str:
    return f'Providing {param!r} is only supported for numpy arrays.'


def is_complex(x: Array, xp: ModuleType) -> bool:
    return xp.isdtype(x.dtype, 'complex floating')


def get_xp_devices(xp: ModuleType) -> list[str] | list[None]:
    """Returns a list of available devices for the given namespace."""
    devices: list[str] = []
    if is_torch(xp):
        devices += ['cpu']
        import torch # type: ignore[import]
        num_cuda = torch.cuda.device_count()
        for i in range(0, num_cuda):
            devices += [f'cuda:{i}']
        if torch.backends.mps.is_available():
            devices += ['mps']
        return devices
    elif is_cupy(xp):
        import cupy # type: ignore[import]
        num_cuda = cupy.cuda.runtime.getDeviceCount()
        for i in range(0, num_cuda):
            devices += [f'cuda:{i}']
        return devices
    elif is_jax(xp):
        import jax # type: ignore[import]
        num_cpu = jax.device_count(backend='cpu')
        for i in range(0, num_cpu):
            devices += [f'cpu:{i}']
        num_gpu = jax.device_count(backend='gpu')
        for i in range(0, num_gpu):
            devices += [f'gpu:{i}']
        num_tpu = jax.device_count(backend='tpu')
        for i in range(0, num_tpu):
            devices += [f'tpu:{i}']
        return devices

    # given namespace is not known to have a list of available devices;
    # return `[None]` so that one can use this in tests for `device=None`.
    return [None]


def scipy_namespace_for(xp: ModuleType) -> ModuleType | None:
    """Return the `scipy`-like namespace of a non-NumPy backend

    That is, return the namespace corresponding with backend `xp` that contains
    `scipy` sub-namespaces like `linalg` and `special`. If no such namespace
    exists, return ``None``. Useful for dispatching.
    """

    if is_cupy(xp):
        import cupyx  # type: ignore[import-not-found,import-untyped]
        return cupyx.scipy

    if is_jax(xp):
        import jax  # type: ignore[import-not-found]
        return jax.scipy

    if is_torch(xp):
        return xp

    return None


# temporary substitute for xp.moveaxis, which is not yet in all backends
# or covered by array_api_compat.
def xp_moveaxis_to_end(
        x: Array,
        source: int,
        /, *,
        xp: ModuleType | None = None) -> Array:
    xp = array_namespace(xp) if xp is None else xp
    axes = list(range(x.ndim))
    temp = axes.pop(source)
    axes = axes + [temp]
    return xp.permute_dims(x, axes)


# temporary substitute for xp.copysign, which is not yet in all backends
# or covered by array_api_compat.
def xp_copysign(x1: Array, x2: Array, /, *, xp: ModuleType | None = None) -> Array:
    # no attempt to account for special cases
    xp = array_namespace(x1, x2) if xp is None else xp
    abs_x1 = xp.abs(x1)
    return xp.where(x2 >= 0, abs_x1, -abs_x1)


# partial substitute for xp.sign, which does not cover the NaN special case
# that I need. (https://github.com/data-apis/array-api-compat/issues/136)
def xp_sign(x: Array, /, *, xp: ModuleType | None = None) -> Array:
    xp = array_namespace(x) if xp is None else xp
    if is_numpy(xp):  # only NumPy implements the special cases correctly
        return xp.sign(x)
    sign = xp.zeros_like(x)
    one = xp.asarray(1, dtype=x.dtype)
    sign = xp.where(x > 0, one, sign)
    sign = xp.where(x < 0, -one, sign)
    sign = xp.where(xp.isnan(x), xp.nan*one, sign)
    return sign

# maybe use `scipy.linalg` if/when array API support is added
def xp_vector_norm(x: Array, /, *,
                   axis: int | tuple[int] | None = None,
                   keepdims: bool = False,
                   ord: int | float = 2,
                   xp: ModuleType | None = None) -> Array:
    xp = array_namespace(x) if xp is None else xp

    if SCIPY_ARRAY_API:
        # check for optional `linalg` extension
        if hasattr(xp, 'linalg'):
            return xp.linalg.vector_norm(x, axis=axis, keepdims=keepdims, ord=ord)
        else:
            if ord != 2:
                raise ValueError(
                    "only the Euclidean norm (`ord=2`) is currently supported in "
                    "`xp_vector_norm` for backends not implementing the `linalg` "
                    "extension."
                )
            # return (x @ x)**0.5
            # or to get the right behavior with nd, complex arrays
            return xp.sum(xp.conj(x) * x, axis=axis, keepdims=keepdims)**0.5
    else:
        # to maintain backwards compatibility
        return np.linalg.norm(x, ord=ord, axis=axis, keepdims=keepdims)


def xp_ravel(x: Array, /, *, xp: ModuleType | None = None) -> Array:
    # Equivalent of np.ravel written in terms of array API
    # Even though it's one line, it comes up so often that it's worth having
    # this function for readability
    xp = array_namespace(x) if xp is None else xp
    return xp.reshape(x, (-1,))


def xp_real(x: Array, /, *, xp: ModuleType | None = None) -> Array:
    # Convenience wrapper of xp.real that allows non-complex input;
    # see data-apis/array-api#824
    xp = array_namespace(x) if xp is None else xp
    return xp.real(x) if xp.isdtype(x.dtype, 'complex floating') else x


def xp_take_along_axis(arr: Array,
                       indices: Array, /, *,
                       axis: int = -1,
                       xp: ModuleType | None = None) -> Array:
    # Dispatcher for np.take_along_axis for backends that support it;
    # see data-apis/array-api/pull#816
    xp = array_namespace(arr) if xp is None else xp
    if is_torch(xp):
        return xp.take_along_dim(arr, indices, dim=axis)
    elif is_array_api_strict(xp):
        raise NotImplementedError("Array API standard does not define take_along_axis")
    else:
        return xp.take_along_axis(arr, indices, axis)


# utility to broadcast arrays and promote to common dtype
def xp_broadcast_promote(*args, ensure_writeable=False, force_floating=False, xp=None):
    xp = array_namespace(*args) if xp is None else xp

    args = [(_asarray(arg, subok=True) if arg is not None else arg) for arg in args]
    args_not_none = [arg for arg in args if arg is not None]

    # determine minimum dtype
    default_float = xp.asarray(1.).dtype
    dtypes = [arg.dtype for arg in args_not_none]
    try:  # follow library's prefered mixed promotion rules
        dtype = xp.result_type(*dtypes)
        if force_floating and xp.isdtype(dtype, 'integral'):
            # If we were to add `default_float` before checking whether the result
            # type is otherwise integral, we risk promotion from lower float.
            dtype = xp.result_type(dtype, default_float)
    except TypeError:  # mixed type promotion isn't defined
        float_dtypes = [dtype for dtype in dtypes
                        if not xp.isdtype(dtype, 'integral')]
        if float_dtypes:
            dtype = xp.result_type(*float_dtypes, default_float)
        elif force_floating:
            dtype = default_float
        else:
            dtype = xp.result_type(*dtypes)

    # determine result shape
    shapes = {arg.shape for arg in args_not_none}
    try:
        shape = (np.broadcast_shapes(*shapes) if len(shapes) != 1
                 else args_not_none[0].shape)
    except ValueError as e:
        message = "Array shapes are incompatible for broadcasting."
        raise ValueError(message) from e

    out = []
    for arg in args:
        if arg is None:
            out.append(arg)
            continue

        # broadcast only if needed
        # Even if two arguments need broadcasting, this is faster than
        # `broadcast_arrays`, especially since we've already determined `shape`
        if arg.shape != shape:
            kwargs = {'subok': True} if is_numpy(xp) else {}
            arg = xp.broadcast_to(arg, shape, **kwargs)

        # convert dtype/copy only if needed
        if (arg.dtype != dtype) or ensure_writeable:
            arg = xp.astype(arg, dtype, copy=True)
        out.append(arg)

    return out


def xp_float_to_complex(arr: Array, xp: ModuleType | None = None) -> Array:
    xp = array_namespace(arr) if xp is None else xp
    arr_dtype = arr.dtype
    # The standard float dtypes are float32 and float64.
    # Convert float32 to complex64,
    # and float64 (and non-standard real dtypes) to complex128
    if xp.isdtype(arr_dtype, xp.float32):
        arr = xp.astype(arr, xp.complex64)
    elif xp.isdtype(arr_dtype, 'real floating'):
        arr = xp.astype(arr, xp.complex128)

    return arr


def xp_default_dtype(xp):
    """Query the namespace-dependent default floating-point dtype.
    """
    if is_torch(xp):
        # historically, we allow pytorch to keep its default of float32
        return xp.get_default_dtype()
    else:
        # we default to float64
        return xp.float64
