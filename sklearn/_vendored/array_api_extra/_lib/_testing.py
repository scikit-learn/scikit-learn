"""
Testing utilities.

Note that this is private API; don't expect it to be stable.
See also ..testing for public testing utilities.
"""

import math
from types import ModuleType
from typing import cast

import pytest

from ._utils._compat import (
    array_namespace,
    is_array_api_strict_namespace,
    is_cupy_namespace,
    is_dask_namespace,
    is_pydata_sparse_namespace,
    is_torch_namespace,
)
from ._utils._typing import Array

__all__ = ["xp_assert_close", "xp_assert_equal"]


def _check_ns_shape_dtype(
    actual: Array, desired: Array
) -> ModuleType:  # numpydoc ignore=RT03
    """
    Assert that namespace, shape and dtype of the two arrays match.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).

    Returns
    -------
    Arrays namespace.
    """
    actual_xp = array_namespace(actual)  # Raises on scalars and lists
    desired_xp = array_namespace(desired)

    msg = f"namespaces do not match: {actual_xp} != f{desired_xp}"
    assert actual_xp == desired_xp, msg

    actual_shape = actual.shape
    desired_shape = desired.shape
    if is_dask_namespace(desired_xp):
        # Dask uses nan instead of None for unknown shapes
        if any(math.isnan(i) for i in cast(tuple[float, ...], actual_shape)):
            actual_shape = actual.compute().shape  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
        if any(math.isnan(i) for i in cast(tuple[float, ...], desired_shape)):
            desired_shape = desired.compute().shape  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]

    msg = f"shapes do not match: {actual_shape} != f{desired_shape}"
    assert actual_shape == desired_shape, msg

    msg = f"dtypes do not match: {actual.dtype} != {desired.dtype}"
    assert actual.dtype == desired.dtype, msg

    return desired_xp


def xp_assert_equal(actual: Array, desired: Array, err_msg: str = "") -> None:
    """
    Array-API compatible version of `np.testing.assert_array_equal`.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    err_msg : str, optional
        Error message to display on failure.

    See Also
    --------
    xp_assert_close : Similar function for inexact equality checks.
    numpy.testing.assert_array_equal : Similar function for NumPy arrays.
    """
    xp = _check_ns_shape_dtype(actual, desired)

    if is_cupy_namespace(xp):
        xp.testing.assert_array_equal(actual, desired, err_msg=err_msg)
    elif is_torch_namespace(xp):
        # PyTorch recommends using `rtol=0, atol=0` like this
        # to test for exact equality
        xp.testing.assert_close(
            actual,
            desired,
            rtol=0,
            atol=0,
            equal_nan=True,
            check_dtype=False,
            msg=err_msg or None,
        )
    else:
        import numpy as np  # pylint: disable=import-outside-toplevel

        if is_pydata_sparse_namespace(xp):
            actual = actual.todense()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
            desired = desired.todense()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]

        actual_np = None
        desired_np = None
        if is_array_api_strict_namespace(xp):
            # __array__ doesn't work on array-api-strict device arrays
            # We need to convert to the CPU device first
            actual_np = np.asarray(xp.asarray(actual, device=xp.Device("CPU_DEVICE")))
            desired_np = np.asarray(xp.asarray(desired, device=xp.Device("CPU_DEVICE")))

        # JAX/Dask arrays work with `np.testing`
        actual_np = actual if actual_np is None else actual_np
        desired_np = desired if desired_np is None else desired_np
        np.testing.assert_array_equal(actual_np, desired_np, err_msg=err_msg)  # pyright: ignore[reportUnknownArgumentType]


def xp_assert_close(
    actual: Array,
    desired: Array,
    *,
    rtol: float | None = None,
    atol: float = 0,
    err_msg: str = "",
) -> None:
    """
    Array-API compatible version of `np.testing.assert_allclose`.

    Parameters
    ----------
    actual : Array
        The array produced by the tested function.
    desired : Array
        The expected array (typically hardcoded).
    rtol : float, optional
        Relative tolerance. Default: dtype-dependent.
    atol : float, optional
        Absolute tolerance. Default: 0.
    err_msg : str, optional
        Error message to display on failure.

    See Also
    --------
    xp_assert_equal : Similar function for exact equality checks.
    isclose : Public function for checking closeness.
    numpy.testing.assert_allclose : Similar function for NumPy arrays.

    Notes
    -----
    The default `atol` and `rtol` differ from `xp.all(xpx.isclose(a, b))`.
    """
    xp = _check_ns_shape_dtype(actual, desired)

    floating = xp.isdtype(actual.dtype, ("real floating", "complex floating"))
    if rtol is None and floating:
        # multiplier of 4 is used as for `np.float64` this puts the default `rtol`
        # roughly half way between sqrt(eps) and the default for
        # `numpy.testing.assert_allclose`, 1e-7
        rtol = xp.finfo(actual.dtype).eps ** 0.5 * 4
    elif rtol is None:
        rtol = 1e-7

    if is_cupy_namespace(xp):
        xp.testing.assert_allclose(
            actual, desired, rtol=rtol, atol=atol, err_msg=err_msg
        )
    elif is_torch_namespace(xp):
        xp.testing.assert_close(
            actual, desired, rtol=rtol, atol=atol, equal_nan=True, msg=err_msg or None
        )
    else:
        import numpy as np  # pylint: disable=import-outside-toplevel

        if is_pydata_sparse_namespace(xp):
            actual = actual.todense()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
            desired = desired.todense()  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]

        actual_np = None
        desired_np = None
        if is_array_api_strict_namespace(xp):
            # __array__ doesn't work on array-api-strict device arrays
            # We need to convert to the CPU device first
            actual_np = np.asarray(xp.asarray(actual, device=xp.Device("CPU_DEVICE")))
            desired_np = np.asarray(xp.asarray(desired, device=xp.Device("CPU_DEVICE")))

        # JAX/Dask arrays work with `np.testing`
        actual_np = actual if actual_np is None else actual_np
        desired_np = desired if desired_np is None else desired_np

        assert isinstance(rtol, float)
        np.testing.assert_allclose(  # pyright: ignore[reportCallIssue]
            actual_np,  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]
            desired_np,  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]
            rtol=rtol,
            atol=atol,
            err_msg=err_msg,
        )


def xfail(request: pytest.FixtureRequest, reason: str) -> None:
    """
    XFAIL the currently running test.

    Unlike ``pytest.xfail``, allow rest of test to execute instead of immediately
    halting it, so that it may result in a XPASS.
    xref https://github.com/pandas-dev/pandas/issues/38902

    Parameters
    ----------
    request : pytest.FixtureRequest
        ``request`` argument of the test function.
    reason : str
        Reason for the expected failure.
    """
    request.node.add_marker(pytest.mark.xfail(reason=reason))
