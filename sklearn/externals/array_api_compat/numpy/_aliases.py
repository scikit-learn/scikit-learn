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
permute_dims = get_xp(np)(_aliases.permute_dims)
if np.__version__ < '2.1':
    reshape = get_xp(np)(_aliases.reshape)
else:
    reshape = np.reshape
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


def clip(
    x: Array,
    /,
    min: float | Array | None = None,
    max: float | Array | None = None,
    **kwargs,
) -> Array:
    """Array API compatible clip implementation for NumPy.

    NumPy's native ``clip`` is used directly after casting bounds to the
    input dtype. This keeps the result dtype aligned with ``x.dtype`` and
    avoids NumPy's default promotion behavior.

    Args:
        x: Input array.
        min: Minimum bound. If None, no lower bound is applied.
        max: Maximum bound. If None, no upper bound is applied.
        **kwargs: Additional keyword arguments passed to ``np.clip``.
        out: Optional output array to store the result, has to have dtype of x
    """
    # out is a possible *kwarg for numpy.clip, but not in the array API spec. We handle it here to
    # avoid having to add it to the array API spec, which would be a breaking change
    # check if out in kwargs, if so pop it and use it as the out parameter
    if "out" in kwargs:
        out = kwargs.pop("out")
    else:
        out = None

    def _bound_shape(a: object) -> tuple[int, ...]:
        if a is None or np.isscalar(a):
            return ()
        return np.asarray(a).shape

    dtype = x.dtype
    out_dtype = out.dtype if out is not None else dtype
    if out_dtype != dtype:
        raise ValueError(
            f"Output array has dtype {out_dtype}, but input array has dtype {dtype}"
            )
    min_shape = _bound_shape(min)
    max_shape = _bound_shape(max)

    # avoid shape broadcasting and copying when not necessary
    if min_shape == () and max_shape == ():
        result_shape = x.shape
    else:
        result_shape = np.broadcast_shapes(x.shape, min_shape, max_shape)

    # Handle cases where the bounds are outside the range of the integer input dtype
    # this covers integer arrays for float and integer bounds
    # Also handle cases where the min/max are arrays/lists (replace values below min with iinfo.min and above max with iinfo.max)
    if np.issubdtype(dtype, np.integer):
        
        if np.issubdtype(type(min), np.integer) and min <= np.iinfo(dtype).min:
            min = None
        elif np.issubdtype(type(min), np.floating) and min < np.iinfo(dtype).min:
            min = np.iinfo(dtype).min
        elif isinstance(min, Array):
            min[min < np.iinfo(dtype).min] = np.iinfo(dtype).min
        elif isinstance(min, (list,tuple)):
            min = np.asarray(min)
            min[min < np.iinfo(dtype).min] = np.iinfo(dtype).min
            
        
        if np.issubdtype(type(max), np.integer) and max >= np.iinfo(dtype).max:
            max = None
        elif np.issubdtype(type(max), np.floating) and max > np.iinfo(dtype).max:
            max = np.iinfo(dtype).max
        elif isinstance(max, Array):
            max[max > np.iinfo(dtype).max] = np.iinfo(dtype).max
        elif isinstance(max, (list,tuple)):
            max = np.asarray(max)
            max[max > np.iinfo(dtype).max] = np.iinfo(dtype).max
    
    # In the case of downcasting floats numpy replaces out of bounds with inf 
    # This automatically handles those cases
    
    # Early return for simple cases
    if min is None and max is None:
        if out is None:
            return x.copy()[()]
        np.copyto(out, x)
        return out[()]

    # Cast clip parameters to the input dtype and broadcast them to the result shape.
    a_min = None
    if min is not None:
        a_min = np.asarray(min, dtype=dtype)
        if a_min.shape != result_shape:
            # Casting first keeps NumPy from promoting the output dtype.
            a_min = np.broadcast_to(a_min, result_shape)

    a_max = None
    if max is not None:
        a_max = np.asarray(max, dtype=dtype)
        if a_max.shape != result_shape:
            # Casting first keeps NumPy from promoting the output dtype.
            a_max = np.broadcast_to(a_max, result_shape)

    if out is None:
        out = np.empty(result_shape, dtype=dtype)

    np.clip(x, a_min, a_max, out=out, casting="no", **kwargs)
    return out[()]


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
    "clip",
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
    "take_along_axis",
]


def __dir__() -> list[str]:
    return __all__
