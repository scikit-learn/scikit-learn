"""Utility functions to use Python Array API compatible libraries.

For the context about the Array API see:
https://data-apis.org/array-api/latest/purpose_and_scope.html

The SciPy use case of the Array API is described on the following page:
https://data-apis.org/array-api/latest/use_cases.html#use-case-scipy
"""
from __future__ import annotations

import os
import warnings

from types import ModuleType
from typing import Any, Literal, TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from scipy._lib import array_api_compat
from scipy._lib.array_api_compat import (
    is_array_api_obj,
    size,
    numpy as np_compat,
    device
)

__all__ = ['array_namespace', '_asarray', 'size', 'device']


# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")

_GLOBAL_CONFIG = {
    "SCIPY_ARRAY_API": SCIPY_ARRAY_API,
    "SCIPY_DEVICE": SCIPY_DEVICE,
}


if TYPE_CHECKING:
    Array = Any  # To be changed to a Protocol later (see array-api#589)
    ArrayLike = Array | npt.ArrayLike


def compliance_scipy(arrays: list[ArrayLike]) -> list[Array]:
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
        if isinstance(array, np.ma.MaskedArray):
            raise TypeError("Inputs of type `numpy.ma.MaskedArray` are not supported.")
        elif isinstance(array, np.matrix):
            raise TypeError("Inputs of type `numpy.matrix` are not supported.")
        if isinstance(array, (np.ndarray, np.generic)):
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
    2. `compliance_scipy` raise exceptions on known-bad subclasses. See
       its definition for more details.

    When the global switch is False, it defaults to the `numpy` namespace.
    In that case, there is no compliance check. This is a convenience to
    ease the adoption. Otherwise, arrays must comply with the new rules.
    """
    if not _GLOBAL_CONFIG["SCIPY_ARRAY_API"]:
        # here we could wrap the namespace if needed
        return np_compat

    _arrays = [array for array in arrays if array is not None]

    _arrays = compliance_scipy(_arrays)

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
    if xp.__name__ in {"numpy", "scipy._lib.array_api_compat.numpy"}:
        # Use NumPy API to support order
        if copy is True:
            array = np.array(array, order=order, dtype=dtype, subok=subok)
        elif subok:
            array = np.asanyarray(array, order=order, dtype=dtype)
        else:
            array = np.asarray(array, order=order, dtype=dtype)

        # At this point array is a NumPy ndarray. We convert it to an array
        # container that is consistent with the input's namespace.
        array = xp.asarray(array)
    else:
        try:
            array = xp.asarray(array, dtype=dtype, copy=copy)
        except TypeError:
            coerced_xp = array_namespace(xp.asarray(3))
            array = coerced_xp.asarray(array, dtype=dtype, copy=copy)

    if check_finite:
        _check_finite(array, xp)

    return array


def atleast_nd(x: Array, *, ndim: int, xp: ModuleType | None = None) -> Array:
    """Recursively expand the dimension to have at least `ndim`."""
    if xp is None:
        xp = array_namespace(x)
    x = xp.asarray(x)
    if x.ndim < ndim:
        x = xp.expand_dims(x, axis=0)
        x = atleast_nd(x, ndim=ndim, xp=xp)
    return x


def copy(x: Array, *, xp: ModuleType | None = None) -> Array:
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
    # Note: xp.asarray fails if xp is numpy.
    if xp is None:
        xp = array_namespace(x)

    return _asarray(x, copy=True, xp=xp)


def is_numpy(xp: ModuleType) -> bool:
    return xp.__name__ in ('numpy', 'scipy._lib.array_api_compat.numpy')


def is_cupy(xp: ModuleType) -> bool:
    return xp.__name__ in ('cupy', 'scipy._lib.array_api_compat.cupy')


def is_torch(xp: ModuleType) -> bool:
    return xp.__name__ in ('torch', 'scipy._lib.array_api_compat.torch')

def is_jax(xp):
    return xp.__name__ in ('jax.numpy', 'jax.experimental.array_api')


def _strict_check(actual, desired, xp,
                  check_namespace=True, check_dtype=True, check_shape=True):
    __tracebackhide__ = True  # Hide traceback for py.test
    if check_namespace:
        _assert_matching_namespace(actual, desired)

    desired = xp.asarray(desired)

    if check_dtype:
        _msg = f"dtypes do not match.\nActual: {actual.dtype}\nDesired: {desired.dtype}"
        assert actual.dtype == desired.dtype, _msg

    if check_shape:
        _msg = f"Shapes do not match.\nActual: {actual.shape}\nDesired: {desired.shape}"
        assert actual.shape == desired.shape, _msg
        _check_scalar(actual, desired, xp)

    desired = xp.broadcast_to(desired, actual.shape)
    return desired


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


def _check_scalar(actual, desired, xp):
    __tracebackhide__ = True  # Hide traceback for py.test
    # Shape check alone is sufficient unless desired.shape == (). Also,
    # only NumPy distinguishes between scalars and arrays.
    if desired.shape != () or not is_numpy(xp):
        return
    # We want to follow the conventions of the `xp` library. Libraries like
    # NumPy, for which `np.asarray(0)[()]` returns a scalar, tend to return
    # a scalar even when a 0D array might be more appropriate:
    # import numpy as np
    # np.mean([1, 2, 3])  # scalar, not 0d array
    # np.asarray(0)*2  # scalar, not 0d array
    # np.sin(np.asarray(0))  # scalar, not 0d array
    # Libraries like CuPy, for which `cp.asarray(0)[()]` returns a 0D array,
    # tend to return a 0D array in scenarios like those above.
    # Therefore, regardless of whether the developer provides a scalar or 0D
    # array for `desired`, we would typically want the type of `actual` to be
    # the type of `desired[()]`. If the developer wants to override this
    # behavior, they can set `check_shape=False`.
    desired = desired[()]
    _msg = f"Types do not match:\n Actual: {type(actual)}\n Desired: {type(desired)}"
    assert (xp.isscalar(actual) and xp.isscalar(desired)
            or (not xp.isscalar(actual) and not xp.isscalar(desired))), _msg


def xp_assert_equal(actual, desired, check_namespace=True, check_dtype=True,
                    check_shape=True, err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)
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


def xp_assert_close(actual, desired, rtol=None, atol=0, check_namespace=True,
                    check_dtype=True, check_shape=True, err_msg='', xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)

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


def xp_assert_less(actual, desired, check_namespace=True, check_dtype=True,
                   check_shape=True, err_msg='', verbose=True, xp=None):
    __tracebackhide__ = True  # Hide traceback for py.test
    if xp is None:
        xp = array_namespace(actual)
    desired = _strict_check(actual, desired, xp, check_namespace=check_namespace,
                            check_dtype=check_dtype, check_shape=check_shape)
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


def cov(x: Array, *, xp: ModuleType | None = None) -> Array:
    if xp is None:
        xp = array_namespace(x)

    X = copy(x, xp=xp)
    dtype = xp.result_type(X, xp.float64)

    X = atleast_nd(X, ndim=2, xp=xp)
    X = xp.asarray(X, dtype=dtype)

    avg = xp.mean(X, axis=1)
    fact = X.shape[1] - 1

    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice",
                      RuntimeWarning, stacklevel=2)
        fact = 0.0

    X -= avg[:, None]
    X_T = X.T
    if xp.isdtype(X_T.dtype, 'complex floating'):
        X_T = xp.conj(X_T)
    c = X @ X_T
    c /= fact
    axes = tuple(axis for axis, length in enumerate(c.shape) if length == 1)
    return xp.squeeze(c, axis=axes)


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


def scipy_namespace_for(xp: ModuleType) -> ModuleType:
    """
    Return the `scipy` namespace for alternative backends, where it exists,
    such as `cupyx.scipy` and `jax.scipy`. Useful for ad hoc dispatching.

    Default: return `scipy` (this package).
    """


    if is_cupy(xp):
        import cupyx  # type: ignore[import-not-found,import-untyped]
        return cupyx.scipy

    if is_jax(xp):
        import jax  # type: ignore[import-not-found]
        return jax.scipy

    import scipy
    return scipy


# temporary substitute for xp.minimum, which is not yet in all backends
# or covered by array_api_compat.
def xp_minimum(x1: Array, x2: Array, /) -> Array:
    # xp won't be passed in because it doesn't need to be passed in to xp.minimum
    xp = array_namespace(x1, x2)
    if hasattr(xp, 'minimum'):
        return xp.minimum(x1, x2)
    x1, x2 = xp.broadcast_arrays(x1, x2)
    i = (x2 < x1) | xp.isnan(x2)
    res = xp.where(i, x2, x1)
    return res[()] if res.ndim == 0 else res


# temporary substitute for xp.clip, which is not yet in all backends
# or covered by array_api_compat.
def xp_clip(
        x: Array,
        /,
        min: int | float | Array | None = None,
        max: int | float | Array | None = None,
        *,
        xp: ModuleType | None = None) -> Array:
    xp = array_namespace(x) if xp is None else xp
    a, b = xp.asarray(min, dtype=x.dtype), xp.asarray(max, dtype=x.dtype)
    if hasattr(xp, 'clip'):
        return xp.clip(x, a, b)
    x, a, b = xp.broadcast_arrays(x, a, b)
    y = xp.asarray(x, copy=True)
    ia = y < a
    y[ia] = a[ia]
    ib = y > b
    y[ib] = b[ib]
    return y[()] if y.ndim == 0 else y


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
    sign = xp.full_like(x, xp.nan)
    one = xp.asarray(1, dtype=x.dtype)
    sign = xp.where(x > 0, one, sign)
    sign = xp.where(x < 0, -one, sign)
    sign = xp.where(x == 0, 0*one, sign)
    return sign
