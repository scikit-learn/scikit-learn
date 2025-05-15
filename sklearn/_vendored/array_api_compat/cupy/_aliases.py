from __future__ import annotations

import cupy as cp

from ..common import _aliases, _helpers
from .._internal import get_xp

from ._info import __array_namespace_info__

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Optional, Union
    from ._typing import ndarray, Device, Dtype, NestedSequence, SupportsBufferProtocol

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
ceil = get_xp(cp)(_aliases.ceil)
floor = get_xp(cp)(_aliases.floor)
trunc = get_xp(cp)(_aliases.trunc)
matmul = get_xp(cp)(_aliases.matmul)
matrix_transpose = get_xp(cp)(_aliases.matrix_transpose)
tensordot = get_xp(cp)(_aliases.tensordot)
sign = get_xp(cp)(_aliases.sign)

_copy_default = object()

# asarray also adds the copy keyword, which is not present in numpy 1.0.
def asarray(
    obj: Union[
        ndarray,
        bool,
        int,
        float,
        NestedSequence[bool | int | float],
        SupportsBufferProtocol,
    ],
    /,
    *,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    copy: Optional[bool] = _copy_default,
    **kwargs,
) -> ndarray:
    """
    Array API compatibility wrapper for asarray().

    See the corresponding documentation in the array library and/or the array API
    specification for more details.
    """
    with cp.cuda.Device(device):
        # cupy is like NumPy 1.26 (except without _CopyMode). See the comments
        # in asarray in numpy/_aliases.py.
        if copy is not _copy_default:
            # A future version of CuPy will change the meaning of copy=False
            # to mean no-copy. We don't know for certain what version it will
            # be yet, so to avoid breaking that version, we use a different
            # default value for copy so asarray(obj) with no copy kwarg will
            # always do the copy-if-needed behavior.

            # This will still need to be updated to remove the
            # NotImplementedError for copy=False, but at least this won't
            # break the default or existing behavior.
            if copy is None:
                copy = False
            elif copy is False:
                raise NotImplementedError("asarray(copy=False) is not yet supported in cupy")
            kwargs['copy'] = copy

        return cp.array(obj, dtype=dtype, **kwargs)


def astype(
    x: ndarray,
    dtype: Dtype,
    /,
    *,
    copy: bool = True,
    device: Optional[Device] = None,
) -> ndarray:
    if device is None:
        return x.astype(dtype=dtype, copy=copy)
    out = _helpers.to_device(x.astype(dtype=dtype, copy=False), device)
    return out.copy() if copy and out is x else out


# cupy.count_nonzero does not have keepdims
def count_nonzero(
    x: ndarray,
    axis=None,
    keepdims=False
) -> ndarray:
   result = cp.count_nonzero(x, axis)
   if keepdims:
       if axis is None:
            return cp.reshape(result, [1]*x.ndim)
       return cp.expand_dims(result, axis)
   return result


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

__all__ = _aliases.__all__ + ['__array_namespace_info__', 'asarray', 'astype',
                              'acos', 'acosh', 'asin', 'asinh', 'atan',
                              'atan2', 'atanh', 'bitwise_left_shift',
                              'bitwise_invert', 'bitwise_right_shift',
                              'bool', 'concat', 'count_nonzero', 'pow', 'sign']

_all_ignore = ['cp', 'get_xp']
