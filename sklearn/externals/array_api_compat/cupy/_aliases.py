from __future__ import annotations

from builtins import bool as py_bool

import cupy as cp

from ..common import _aliases, _helpers
from ..common._typing import NestedSequence, SupportsBufferProtocol
from .._internal import get_xp
from ._typing import Array, Device, DType

bool = cp.bool_

# Basic renames
acos = cp.arccos
acosh = cp.arccosh
asin = cp.arcsin
asinh = cp.arcsinh
atan = cp.arctan
atan2 = cp.arctan2
atanh = cp.arctanh
bitwise_left_shift = cp.left_shift
bitwise_invert = cp.invert
bitwise_right_shift = cp.right_shift
concat = cp.concatenate
pow = cp.power

arange = get_xp(cp)(_aliases.arange)
empty = get_xp(cp)(_aliases.empty)
empty_like = get_xp(cp)(_aliases.empty_like)
eye = get_xp(cp)(_aliases.eye)
full = get_xp(cp)(_aliases.full)
full_like = get_xp(cp)(_aliases.full_like)
linspace = get_xp(cp)(_aliases.linspace)
ones = get_xp(cp)(_aliases.ones)
ones_like = get_xp(cp)(_aliases.ones_like)
zeros = get_xp(cp)(_aliases.zeros)
zeros_like = get_xp(cp)(_aliases.zeros_like)
UniqueAllResult = get_xp(cp)(_aliases.UniqueAllResult)
UniqueCountsResult = get_xp(cp)(_aliases.UniqueCountsResult)
UniqueInverseResult = get_xp(cp)(_aliases.UniqueInverseResult)
unique_all = get_xp(cp)(_aliases.unique_all)
unique_counts = get_xp(cp)(_aliases.unique_counts)
unique_inverse = get_xp(cp)(_aliases.unique_inverse)
unique_values = get_xp(cp)(_aliases.unique_values)
std = get_xp(cp)(_aliases.std)
var = get_xp(cp)(_aliases.var)
cumulative_sum = get_xp(cp)(_aliases.cumulative_sum)
cumulative_prod = get_xp(cp)(_aliases.cumulative_prod)
clip = get_xp(cp)(_aliases.clip)
permute_dims = get_xp(cp)(_aliases.permute_dims)
reshape = get_xp(cp)(_aliases.reshape)
argsort = get_xp(cp)(_aliases.argsort)
sort = get_xp(cp)(_aliases.sort)
nonzero = get_xp(cp)(_aliases.nonzero)
matmul = get_xp(cp)(_aliases.matmul)
matrix_transpose = get_xp(cp)(_aliases.matrix_transpose)
tensordot = get_xp(cp)(_aliases.tensordot)
sign = get_xp(cp)(_aliases.sign)
finfo = get_xp(cp)(_aliases.finfo)
iinfo = get_xp(cp)(_aliases.iinfo)


# asarray also adds the copy keyword, which is not present in numpy 1.0.
def asarray(
    obj: Array | complex | NestedSequence[complex] | SupportsBufferProtocol,
    /,
    *,
    dtype: DType | None = None,
    device: Device | None = None,
    copy: py_bool | None = None,
    **kwargs: object,
) -> Array:
    """
    Array API compatibility wrapper for asarray().

    See the corresponding documentation in the array library and/or the array API
    specification for more details.
    """
    with cp.cuda.Device(device):
        if copy is None:
            return cp.asarray(obj, dtype=dtype, **kwargs)
        else:
            res = cp.array(obj, dtype=dtype, copy=copy, **kwargs)
            if not copy and res is not obj:
                raise ValueError("Unable to avoid copy while creating an array as requested")
            return res


def astype(
    x: Array,
    dtype: DType,
    /,
    *,
    copy: py_bool = True,
    device: Device | None = None,
) -> Array:
    if device is None:
        return x.astype(dtype=dtype, copy=copy)
    out = _helpers.to_device(x.astype(dtype=dtype, copy=False), device)
    return out.copy() if copy and out is x else out


# cupy.count_nonzero does not have keepdims
def count_nonzero(
    x: Array,
    axis: int | tuple[int, ...] | None = None,
    keepdims: py_bool = False,
) -> Array:
   result = cp.count_nonzero(x, axis)
   if keepdims:
       if axis is None:
            return cp.reshape(result, [1]*x.ndim)
       return cp.expand_dims(result, axis)
   return result

# ceil, floor, and trunc return integers for integer inputs

def ceil(x: Array, /) -> Array:
    if cp.issubdtype(x.dtype, cp.integer):
        return x.copy()
    return cp.ceil(x)


def floor(x: Array, /) -> Array:
    if cp.issubdtype(x.dtype, cp.integer):
        return x.copy()
    return cp.floor(x)


def trunc(x: Array, /) -> Array:
    if cp.issubdtype(x.dtype, cp.integer):
        return x.copy()
    return cp.trunc(x)


# take_along_axis: axis defaults to -1 but in cupy (and numpy) axis is a required arg
def take_along_axis(x: Array, indices: Array, /, *, axis: int = -1) -> Array:
    return cp.take_along_axis(x, indices, axis=axis)


# These functions are completely new here. If the library already has them
# (i.e., numpy 2.0), use the library version instead of our wrapper.
if hasattr(cp, 'vecdot'):
    vecdot = cp.vecdot
else:
    vecdot = get_xp(cp)(_aliases.vecdot)

if hasattr(cp, 'isdtype'):
    isdtype = cp.isdtype
else:
    isdtype = get_xp(cp)(_aliases.isdtype)

if hasattr(cp, 'unstack'):
    unstack = cp.unstack
else:
    unstack = get_xp(cp)(_aliases.unstack)

__all__ = _aliases.__all__ + ['asarray', 'astype',
                              'acos', 'acosh', 'asin', 'asinh', 'atan',
                              'atan2', 'atanh', 'bitwise_left_shift',
                              'bitwise_invert', 'bitwise_right_shift',
                              'bool', 'concat', 'count_nonzero', 'pow', 'sign',
                              'ceil', 'floor', 'trunc', 'take_along_axis']


def __dir__() -> list[str]:
    return __all__
