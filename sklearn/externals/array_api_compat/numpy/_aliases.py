# pyright: reportPrivateUsage=false
from __future__ import annotations

from builtins import bool as py_bool
from typing import Any, cast

import numpy as np

from .._internal import get_xp
from ..common import _aliases, _helpers
from ..common._typing import NestedSequence, SupportsBufferProtocol
from ._typing import Array, Device, DType

bool = np.bool_

# Basic renames
acos = np.arccos
acosh = np.arccosh
asin = np.arcsin
asinh = np.arcsinh
atan = np.arctan
atan2 = np.arctan2
atanh = np.arctanh
bitwise_left_shift = np.left_shift
bitwise_invert = np.invert
bitwise_right_shift = np.right_shift
concat = np.concatenate
pow = np.power

arange = get_xp(np)(_aliases.arange)
empty = get_xp(np)(_aliases.empty)
empty_like = get_xp(np)(_aliases.empty_like)
eye = get_xp(np)(_aliases.eye)
full = get_xp(np)(_aliases.full)
full_like = get_xp(np)(_aliases.full_like)
linspace = get_xp(np)(_aliases.linspace)
ones = get_xp(np)(_aliases.ones)
ones_like = get_xp(np)(_aliases.ones_like)
zeros = get_xp(np)(_aliases.zeros)
zeros_like = get_xp(np)(_aliases.zeros_like)
UniqueAllResult = get_xp(np)(_aliases.UniqueAllResult)
UniqueCountsResult = get_xp(np)(_aliases.UniqueCountsResult)
UniqueInverseResult = get_xp(np)(_aliases.UniqueInverseResult)
unique_all = get_xp(np)(_aliases.unique_all)
unique_counts = get_xp(np)(_aliases.unique_counts)
unique_inverse = get_xp(np)(_aliases.unique_inverse)
unique_values = get_xp(np)(_aliases.unique_values)
std = get_xp(np)(_aliases.std)
var = get_xp(np)(_aliases.var)
cumulative_sum = get_xp(np)(_aliases.cumulative_sum)
cumulative_prod = get_xp(np)(_aliases.cumulative_prod)
clip = get_xp(np)(_aliases.clip)
permute_dims = get_xp(np)(_aliases.permute_dims)
reshape = get_xp(np)(_aliases.reshape)
argsort = get_xp(np)(_aliases.argsort)
sort = get_xp(np)(_aliases.sort)
nonzero = get_xp(np)(_aliases.nonzero)
matmul = get_xp(np)(_aliases.matmul)
matrix_transpose = get_xp(np)(_aliases.matrix_transpose)
tensordot = get_xp(np)(_aliases.tensordot)
sign = get_xp(np)(_aliases.sign)
finfo = get_xp(np)(_aliases.finfo)
iinfo = get_xp(np)(_aliases.iinfo)


# asarray also adds the copy keyword, which is not present in numpy 1.0.
# asarray() is different enough between numpy, cupy, and dask, the logic
# complicated enough that it's easier to define it separately for each module
# rather than trying to combine everything into one function in common/
def asarray(
    obj: Array | complex | NestedSequence[complex] | SupportsBufferProtocol,
    /,
    *,
    dtype: DType | None = None,
    device: Device | None = None,
    copy: py_bool | None = None,
    **kwargs: Any,
) -> Array:
    """
    Array API compatibility wrapper for asarray().

    See the corresponding documentation in the array library and/or the array API
    specification for more details.
    """
    _helpers._check_device(np, device)

    # None is unsupported in NumPy 1.0, but we can use an internal enum
    # False in NumPy 1.0 means None in NumPy 2.0 and in the Array API
    if copy is None:
        copy = np._CopyMode.IF_NEEDED  # type: ignore[assignment,attr-defined]
    elif copy is False:
        copy = np._CopyMode.NEVER  # type: ignore[assignment,attr-defined]

    return np.array(obj, copy=copy, dtype=dtype, **kwargs)


def astype(
    x: Array,
    dtype: DType,
    /,
    *,
    copy: py_bool = True,
    device: Device | None = None,
) -> Array:
    _helpers._check_device(np, device)
    return x.astype(dtype=dtype, copy=copy)


# count_nonzero returns a python int for axis=None and keepdims=False
# https://github.com/numpy/numpy/issues/17562
def count_nonzero(
    x: Array,
    axis: int | tuple[int, ...] | None = None,
    keepdims: py_bool = False,
) -> Array:
    # NOTE: this is currently incorrectly typed in numpy, but will be fixed in
    # numpy 2.2.5 and 2.3.0: https://github.com/numpy/numpy/pull/28750
    result = cast("Any", np.count_nonzero(x, axis=axis, keepdims=keepdims))  # pyright: ignore[reportArgumentType, reportCallIssue]
    if axis is None and not keepdims:
        return np.asarray(result)
    return result


# take_along_axis: axis defaults to -1 but in numpy axis is a required arg
def take_along_axis(x: Array, indices: Array, /, *, axis: int = -1) -> Array:
    return np.take_along_axis(x, indices, axis=axis)


# ceil, floor, and trunc return integers for integer inputs in NumPy < 2

def ceil(x: Array, /) -> Array:
    if np.__version__ < '2' and np.issubdtype(x.dtype, np.integer):
        return x.copy()
    return np.ceil(x)


def floor(x: Array, /) -> Array:
    if np.__version__ < '2' and np.issubdtype(x.dtype, np.integer):
        return x.copy()
    return np.floor(x)


def trunc(x: Array, /) -> Array:
    if np.__version__ < '2' and np.issubdtype(x.dtype, np.integer):
        return x.copy()
    return np.trunc(x)


# These functions are completely new here. If the library already has them
# (i.e., numpy 2.0), use the library version instead of our wrapper.
if hasattr(np, "vecdot"):
    vecdot = np.vecdot
else:
    vecdot = get_xp(np)(_aliases.vecdot)  # type: ignore[assignment]

if hasattr(np, "isdtype"):
    isdtype = np.isdtype
else:
    isdtype = get_xp(np)(_aliases.isdtype)

if hasattr(np, "unstack"):
    unstack = np.unstack
else:
    unstack = get_xp(np)(_aliases.unstack)

__all__ = _aliases.__all__ + [
    "asarray",
    "astype",
    "acos",
    "acosh",
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "ceil",
    "floor",
    "trunc",
    "bitwise_left_shift",
    "bitwise_invert",
    "bitwise_right_shift",
    "bool",
    "concat",
    "count_nonzero",
    "pow",
    "take_along_axis"
]


def __dir__() -> list[str]:
    return __all__
