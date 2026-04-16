"""
Override functions from array_api_compat, for use by array-api-extra
and internally.

See also _array_api_compat_vendor.py
"""
import enum
import os

from functools import lru_cache
from types import ModuleType
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt

from scipy._lib import array_api_compat
import scipy._lib.array_api_compat.numpy as np_compat
from scipy._lib.array_api_compat import is_array_api_obj, is_jax_array
from scipy._lib._sparse import SparseABC


Array: TypeAlias = Any  # To be changed to a Protocol later (see array-api#589)
ArrayLike: TypeAlias = Array | npt.ArrayLike

# To enable array API and strict array-like input validation
SCIPY_ARRAY_API: str | bool = os.environ.get("SCIPY_ARRAY_API", False)
# To control the default device - for use in the test suite only
SCIPY_DEVICE = os.environ.get("SCIPY_DEVICE", "cpu")


class _ArrayClsInfo(enum.Enum):
    skip = 0
    numpy = 1
    array_like = 2
    unknown = 3


@lru_cache(100)
def _validate_array_cls(cls: type) -> _ArrayClsInfo:
    if issubclass(cls, (list,  tuple)):
        return _ArrayClsInfo.array_like

    # this comes from `_util._asarray_validated`
    if issubclass(cls, SparseABC):
        msg = ('Sparse arrays/matrices are not supported by this function. '
                'Perhaps one of the `scipy.sparse.linalg` functions '
                'would work instead.')
        raise ValueError(msg)

    if issubclass(cls, np.ma.MaskedArray):
        raise TypeError("Inputs of type `numpy.ma.MaskedArray` are not supported.")

    if issubclass(cls, np.matrix):
        raise TypeError("Inputs of type `numpy.matrix` are not supported.")

    if issubclass(cls, (np.ndarray, np.generic)):
        return _ArrayClsInfo.numpy

    # Note: this must happen after the test for np.generic, because
    # np.float64 and np.complex128 are subclasses of float and complex respectively.
    # This matches the behavior of array_api_compat.
    if issubclass(cls, (int, float, complex, bool, type(None))):
        return _ArrayClsInfo.skip

    return _ArrayClsInfo.unknown


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
    Wrapper around `array_api_compat.array_namespace`.

    1. Check for the global switch `SCIPY_ARRAY_API`. If disabled, just
       return array_api_compat.numpy namespace and skip all compliance checks.

    2. Check for known-bad array classes.
       The following subclasses are not supported and raise and error:

       - `numpy.ma.MaskedArray`
       - `numpy.matrix`
       - NumPy arrays which do not have a boolean or numerical dtype
       - `scipy.sparse` arrays

    3. Coerce array-likes to NumPy arrays and check their dtype.
       Note that non-scalar array-likes can't be mixed with non-NumPy Array
       API objects; e.g.

       - `array_namespace([1, 2])` returns NumPy namespace;
       - `array_namespace(np.asarray([1, 2], [3, 4])` returns NumPy namespace;
       - `array_namespace(cp.asarray([1, 2], [3, 4])` raises an error.
    """
    if not SCIPY_ARRAY_API:
        # here we could wrap the namespace if needed
        return np_compat

    numpy_arrays = []
    api_arrays = []

    for array in arrays:
        arr_info = _validate_array_cls(type(array))
        if arr_info is _ArrayClsInfo.skip:
            pass

        elif arr_info is _ArrayClsInfo.numpy:
            if array.dtype.kind in 'iufcb':  # Numeric or bool
                numpy_arrays.append(array)
            elif array.dtype.kind == 'V' and is_jax_array(array):
                # Special case for JAX zero gradient arrays;
                # see array_api_compat._common._helpers._is_jax_zero_gradient_array
                api_arrays.append(array)  # JAX zero gradient array
            else:
                raise TypeError(f"An argument has dtype `{array.dtype!r}`; "
                                "only boolean and numerical dtypes are supported.")

        elif arr_info is _ArrayClsInfo.unknown and is_array_api_obj(array):
            api_arrays.append(array)

        else:
            # list, tuple, or arbitrary object
            try:
                array = np.asanyarray(array)
            except TypeError:
                raise TypeError("An argument is neither array API compatible nor "
                                "coercible by NumPy.")
            if array.dtype.kind not in 'iufcb':  # Numeric or bool
                raise TypeError(f"An argument has dtype `{array.dtype!r}`; "
                                "only boolean and numerical dtypes are supported.")
            numpy_arrays.append(array)

    # When there are exclusively NumPy and ArrayLikes, skip calling
    # array_api_compat.array_namespace for performance.
    if not api_arrays:
        return np_compat

    # In case of mix of NumPy/ArrayLike and non-NumPy Array API arrays,
    # let array_api_compat.array_namespace raise an error.
    return array_api_compat.array_namespace(*numpy_arrays, *api_arrays)
