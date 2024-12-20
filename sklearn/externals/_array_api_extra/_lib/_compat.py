### Helpers borrowed from array-api-compat

from __future__ import annotations  # https://github.com/pylint-dev/pylint/pull/9990

import inspect
import sys
import typing

from ._typing import override

if typing.TYPE_CHECKING:
    from ._typing import Array, Device

__all__ = ["device"]


# Placeholder object to represent the dask device
# when the array backend is not the CPU.
# (since it is not easy to tell which device a dask array is on)
class _dask_device:  # pylint: disable=invalid-name
    @override
    def __repr__(self) -> str:
        return "DASK_DEVICE"


_DASK_DEVICE = _dask_device()


# device() is not on numpy.ndarray or dask.array and to_device() is not on numpy.ndarray
# or cupy.ndarray. They are not included in array objects of this library
# because this library just reuses the respective ndarray classes without
# wrapping or subclassing them. These helper functions can be used instead of
# the wrapper functions for libraries that need to support both NumPy/CuPy and
# other libraries that use devices.
def device(x: Array, /) -> Device:
    """
    Hardware device the array data resides on.

    This is equivalent to `x.device` according to the `standard
    <https://data-apis.org/array-api/latest/API_specification/generated/array_api.array.device.html>`__.
    This helper is included because some array libraries either do not have
    the `device` attribute or include it with an incompatible API.

    Parameters
    ----------
    x: array
        array instance from an array API compatible library.

    Returns
    -------
    out: device
        a ``device`` object (see the `Device Support <https://data-apis.org/array-api/latest/design_topics/device_support.html>`__
        section of the array API specification).

    Notes
    -----

    For NumPy the device is always `"cpu"`. For Dask, the device is always a
    special `DASK_DEVICE` object.

    See Also
    --------

    to_device : Move array data to a different device.

    """
    if _is_numpy_array(x):
        return "cpu"
    if _is_dask_array(x):
        # Peek at the metadata of the jax array to determine type
        try:
            import numpy as np  # pylint: disable=import-outside-toplevel

            if isinstance(x._meta, np.ndarray):  # pylint: disable=protected-access
                # Must be on CPU since backed by numpy
                return "cpu"
        except ImportError:
            pass
        return _DASK_DEVICE
    if _is_jax_array(x):
        # JAX has .device() as a method, but it is being deprecated so that it
        # can become a property, in accordance with the standard. In order for
        # this function to not break when JAX makes the flip, we check for
        # both here.
        if inspect.ismethod(x.device):
            return x.device()
        return x.device
    if _is_pydata_sparse_array(x):
        # `sparse` will gain `.device`, so check for this first.
        x_device = getattr(x, "device", None)
        if x_device is not None:
            return x_device
        # Everything but DOK has this attr.
        try:
            inner = x.data
        except AttributeError:
            return "cpu"
        # Return the device of the constituent array
        return device(inner)
    return x.device


def _is_numpy_array(x: Array) -> bool:
    """Return True if `x` is a NumPy array."""
    # Avoid importing NumPy if it isn't already
    if "numpy" not in sys.modules:
        return False

    import numpy as np  # pylint: disable=import-outside-toplevel

    # TODO: Should we reject ndarray subclasses?
    return isinstance(x, (np.ndarray, np.generic)) and not _is_jax_zero_gradient_array(
        x
    )


def _is_dask_array(x: Array) -> bool:
    """Return True if `x` is a dask.array Array."""
    # Avoid importing dask if it isn't already
    if "dask.array" not in sys.modules:
        return False

    # pylint: disable=import-error, import-outside-toplevel
    import dask.array  # type: ignore[import-not-found]  # pyright: ignore[reportMissingImports]

    return isinstance(x, dask.array.Array)


def _is_jax_zero_gradient_array(x: Array) -> bool:
    """Return True if `x` is a zero-gradient array.

    These arrays are a design quirk of Jax that may one day be removed.
    See https://github.com/google/jax/issues/20620.
    """
    if "numpy" not in sys.modules or "jax" not in sys.modules:
        return False

    # pylint: disable=import-error, import-outside-toplevel
    import jax  # type: ignore[import-not-found]  # pyright: ignore[reportMissingImports]
    import numpy as np  # pylint: disable=import-outside-toplevel

    return isinstance(x, np.ndarray) and x.dtype == jax.float0  # pyright: ignore[reportUnknownVariableType]


def _is_jax_array(x: Array) -> bool:
    """Return True if `x` is a JAX array."""
    # Avoid importing jax if it isn't already
    if "jax" not in sys.modules:
        return False

    # pylint: disable=import-error, import-outside-toplevel
    import jax  # pyright: ignore[reportMissingImports]

    return isinstance(x, jax.Array) or _is_jax_zero_gradient_array(x)


def _is_pydata_sparse_array(x: Array) -> bool:
    """Return True if `x` is an array from the `sparse` package."""

    # Avoid importing jax if it isn't already
    if "sparse" not in sys.modules:
        return False

    # pylint: disable=import-error, import-outside-toplevel
    import sparse  # type: ignore[import-not-found]  # pyright: ignore[reportMissingImports]

    # TODO: Account for other backends.
    return isinstance(x, sparse.SparseArray)
