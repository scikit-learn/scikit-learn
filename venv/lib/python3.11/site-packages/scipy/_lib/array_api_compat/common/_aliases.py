"""
These are functions that are just aliases of existing functions in NumPy.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Optional, Sequence, Tuple, Union
    from ._typing import ndarray, Device, Dtype

from typing import NamedTuple
import inspect

from ._helpers import array_namespace, _check_device, device, is_torch_array, is_cupy_namespace

# These functions are modified from the NumPy versions.

# Creation functions add the device keyword (which does nothing for NumPy)

def arange(
    start: Union[int, float],
    /,
    stop: Optional[Union[int, float]] = None,
    step: Union[int, float] = 1,
    *,
    xp,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs
) -> ndarray:
    _check_device(xp, device)
    return xp.arange(start, stop=stop, step=step, dtype=dtype, **kwargs)

def empty(
    shape: Union[int, Tuple[int, ...]],
    xp,
    *,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs
) -> ndarray:
    _check_device(xp, device)
    return xp.empty(shape, dtype=dtype, **kwargs)

def empty_like(
    x: ndarray, /, xp, *, dtype: Optional[Dtype] = None, device: Optional[Device] = None,
    **kwargs
) -> ndarray:
    _check_device(xp, device)
    return xp.empty_like(x, dtype=dtype, **kwargs)

def eye(
    n_rows: int,
    n_cols: Optional[int] = None,
    /,
    *,
    xp,
    k: int = 0,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.eye(n_rows, M=n_cols, k=k, dtype=dtype, **kwargs)

def full(
    shape: Union[int, Tuple[int, ...]],
    fill_value: Union[int, float],
    xp,
    *,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.full(shape, fill_value, dtype=dtype, **kwargs)

def full_like(
    x: ndarray,
    /,
    fill_value: Union[int, float],
    *,
    xp,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.full_like(x, fill_value, dtype=dtype, **kwargs)

def linspace(
    start: Union[int, float],
    stop: Union[int, float],
    /,
    num: int,
    *,
    xp,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    endpoint: bool = True,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.linspace(start, stop, num, dtype=dtype, endpoint=endpoint, **kwargs)

def ones(
    shape: Union[int, Tuple[int, ...]],
    xp,
    *,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.ones(shape, dtype=dtype, **kwargs)

def ones_like(
    x: ndarray, /, xp, *, dtype: Optional[Dtype] = None, device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.ones_like(x, dtype=dtype, **kwargs)

def zeros(
    shape: Union[int, Tuple[int, ...]],
    xp,
    *,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.zeros(shape, dtype=dtype, **kwargs)

def zeros_like(
    x: ndarray, /, xp, *, dtype: Optional[Dtype] = None, device: Optional[Device] = None,
    **kwargs,
) -> ndarray:
    _check_device(xp, device)
    return xp.zeros_like(x, dtype=dtype, **kwargs)

# np.unique() is split into four functions in the array API:
# unique_all, unique_counts, unique_inverse, and unique_values (this is done
# to remove polymorphic return types).

# The functions here return namedtuples (np.unique() returns a normal
# tuple).

# Note that these named tuples aren't actually part of the standard namespace,
# but I don't see any issue with exporting the names here regardless.
class UniqueAllResult(NamedTuple):
    values: ndarray
    indices: ndarray
    inverse_indices: ndarray
    counts: ndarray


class UniqueCountsResult(NamedTuple):
    values: ndarray
    counts: ndarray


class UniqueInverseResult(NamedTuple):
    values: ndarray
    inverse_indices: ndarray


def _unique_kwargs(xp):
    # Older versions of NumPy and CuPy do not have equal_nan. Rather than
    # trying to parse version numbers, just check if equal_nan is in the
    # signature.
    s = inspect.signature(xp.unique)
    if 'equal_nan' in s.parameters:
        return {'equal_nan': False}
    return {}

def unique_all(x: ndarray, /, xp) -> UniqueAllResult:
    kwargs = _unique_kwargs(xp)
    values, indices, inverse_indices, counts = xp.unique(
        x,
        return_counts=True,
        return_index=True,
        return_inverse=True,
        **kwargs,
    )
    # np.unique() flattens inverse indices, but they need to share x's shape
    # See https://github.com/numpy/numpy/issues/20638
    inverse_indices = inverse_indices.reshape(x.shape)
    return UniqueAllResult(
        values,
        indices,
        inverse_indices,
        counts,
    )


def unique_counts(x: ndarray, /, xp) -> UniqueCountsResult:
    kwargs = _unique_kwargs(xp)
    res = xp.unique(
        x,
        return_counts=True,
        return_index=False,
        return_inverse=False,
        **kwargs
    )

    return UniqueCountsResult(*res)


def unique_inverse(x: ndarray, /, xp) -> UniqueInverseResult:
    kwargs = _unique_kwargs(xp)
    values, inverse_indices = xp.unique(
        x,
        return_counts=False,
        return_index=False,
        return_inverse=True,
        **kwargs,
    )
    # xp.unique() flattens inverse indices, but they need to share x's shape
    # See https://github.com/numpy/numpy/issues/20638
    inverse_indices = inverse_indices.reshape(x.shape)
    return UniqueInverseResult(values, inverse_indices)


def unique_values(x: ndarray, /, xp) -> ndarray:
    kwargs = _unique_kwargs(xp)
    return xp.unique(
        x,
        return_counts=False,
        return_index=False,
        return_inverse=False,
        **kwargs,
    )

def astype(x: ndarray, dtype: Dtype, /, *, copy: bool = True) -> ndarray:
    if not copy and dtype == x.dtype:
        return x
    return x.astype(dtype=dtype, copy=copy)

# These functions have different keyword argument names

def std(
    x: ndarray,
    /,
    xp,
    *,
    axis: Optional[Union[int, Tuple[int, ...]]] = None,
    correction: Union[int, float] = 0.0, # correction instead of ddof
    keepdims: bool = False,
    **kwargs,
) -> ndarray:
    return xp.std(x, axis=axis, ddof=correction, keepdims=keepdims, **kwargs)

def var(
    x: ndarray,
    /,
    xp,
    *,
    axis: Optional[Union[int, Tuple[int, ...]]] = None,
    correction: Union[int, float] = 0.0, # correction instead of ddof
    keepdims: bool = False,
    **kwargs,
) -> ndarray:
    return xp.var(x, axis=axis, ddof=correction, keepdims=keepdims, **kwargs)

# cumulative_sum is renamed from cumsum, and adds the include_initial keyword
# argument

def cumulative_sum(
    x: ndarray,
    /,
    xp,
    *,
    axis: Optional[int] = None,
    dtype: Optional[Dtype] = None,
    include_initial: bool = False,
    **kwargs
) -> ndarray:
    wrapped_xp = array_namespace(x)

    # TODO: The standard is not clear about what should happen when x.ndim == 0.
    if axis is None:
        if x.ndim > 1:
            raise ValueError("axis must be specified in cumulative_sum for more than one dimension")
        axis = 0

    res = xp.cumsum(x, axis=axis, dtype=dtype, **kwargs)

    # np.cumsum does not support include_initial
    if include_initial:
        initial_shape = list(x.shape)
        initial_shape[axis] = 1
        res = xp.concatenate(
            [wrapped_xp.zeros(shape=initial_shape, dtype=res.dtype, device=device(res)), res],
            axis=axis,
        )
    return res

# The min and max argument names in clip are different and not optional in numpy, and type
# promotion behavior is different.
def clip(
    x: ndarray,
    /,
    min: Optional[Union[int, float, ndarray]] = None,
    max: Optional[Union[int, float, ndarray]] = None,
    *,
    xp,
    # TODO: np.clip has other ufunc kwargs
    out: Optional[ndarray] = None,
) -> ndarray:
    def _isscalar(a):
        return isinstance(a, (int, float, type(None)))
    min_shape = () if _isscalar(min) else min.shape
    max_shape = () if _isscalar(max) else max.shape

    wrapped_xp = array_namespace(x)

    result_shape = xp.broadcast_shapes(x.shape, min_shape, max_shape)

    # np.clip does type promotion but the array API clip requires that the
    # output have the same dtype as x. We do this instead of just downcasting
    # the result of xp.clip() to handle some corner cases better (e.g.,
    # avoiding uint64 -> float64 promotion).

    # Note: cases where min or max overflow (integer) or round (float) in the
    # wrong direction when downcasting to x.dtype are unspecified. This code
    # just does whatever NumPy does when it downcasts in the assignment, but
    # other behavior could be preferred, especially for integers. For example,
    # this code produces:

    # >>> clip(asarray(0, dtype=int8), asarray(128, dtype=int16), None)
    # -128

    # but an answer of 0 might be preferred. See
    # https://github.com/numpy/numpy/issues/24976 for more discussion on this issue.


    # At least handle the case of Python integers correctly (see
    # https://github.com/numpy/numpy/pull/26892).
    if type(min) is int and min <= wrapped_xp.iinfo(x.dtype).min:
        min = None
    if type(max) is int and max >= wrapped_xp.iinfo(x.dtype).max:
        max = None

    if out is None:
        out = wrapped_xp.asarray(xp.broadcast_to(x, result_shape),
                                 copy=True, device=device(x))
    if min is not None:
        if is_torch_array(x) and x.dtype == xp.float64 and _isscalar(min):
            # Avoid loss of precision due to torch defaulting to float32
            min = wrapped_xp.asarray(min, dtype=xp.float64)
        a = xp.broadcast_to(wrapped_xp.asarray(min, device=device(x)), result_shape)
        ia = (out < a) | xp.isnan(a)
        # torch requires an explicit cast here
        out[ia] = wrapped_xp.astype(a[ia], out.dtype)
    if max is not None:
        if is_torch_array(x) and x.dtype == xp.float64 and _isscalar(max):
            max = wrapped_xp.asarray(max, dtype=xp.float64)
        b = xp.broadcast_to(wrapped_xp.asarray(max, device=device(x)), result_shape)
        ib = (out > b) | xp.isnan(b)
        out[ib] = wrapped_xp.astype(b[ib], out.dtype)
    # Return a scalar for 0-D
    return out[()]

# Unlike transpose(), the axes argument to permute_dims() is required.
def permute_dims(x: ndarray, /, axes: Tuple[int, ...], xp) -> ndarray:
    return xp.transpose(x, axes)

# np.reshape calls the keyword argument 'newshape' instead of 'shape'
def reshape(x: ndarray,
            /,
            shape: Tuple[int, ...],
            xp, copy: Optional[bool] = None,
            **kwargs) -> ndarray:
    if copy is True:
        x = x.copy()
    elif copy is False:
        y = x.view()
        y.shape = shape
        return y
    return xp.reshape(x, shape, **kwargs)

# The descending keyword is new in sort and argsort, and 'kind' replaced with
# 'stable'
def argsort(
    x: ndarray, /, xp, *, axis: int = -1, descending: bool = False, stable: bool = True,
    **kwargs,
) -> ndarray:
    # Note: this keyword argument is different, and the default is different.
    # We set it in kwargs like this because numpy.sort uses kind='quicksort'
    # as the default whereas cupy.sort uses kind=None.
    if stable:
        kwargs['kind'] = "stable"
    if not descending:
        res = xp.argsort(x, axis=axis, **kwargs)
    else:
        # As NumPy has no native descending sort, we imitate it here. Note that
        # simply flipping the results of xp.argsort(x, ...) would not
        # respect the relative order like it would in native descending sorts.
        res = xp.flip(
            xp.argsort(xp.flip(x, axis=axis), axis=axis, **kwargs),
            axis=axis,
        )
        # Rely on flip()/argsort() to validate axis
        normalised_axis = axis if axis >= 0 else x.ndim + axis
        max_i = x.shape[normalised_axis] - 1
        res = max_i - res
    return res

def sort(
    x: ndarray, /, xp, *, axis: int = -1, descending: bool = False, stable: bool = True,
    **kwargs,
) -> ndarray:
    # Note: this keyword argument is different, and the default is different.
    # We set it in kwargs like this because numpy.sort uses kind='quicksort'
    # as the default whereas cupy.sort uses kind=None.
    if stable:
        kwargs['kind'] = "stable"
    res = xp.sort(x, axis=axis, **kwargs)
    if descending:
        res = xp.flip(res, axis=axis)
    return res

# nonzero should error for zero-dimensional arrays
def nonzero(x: ndarray, /, xp, **kwargs) -> Tuple[ndarray, ...]:
    if x.ndim == 0:
        raise ValueError("nonzero() does not support zero-dimensional arrays")
    return xp.nonzero(x, **kwargs)

# ceil, floor, and trunc return integers for integer inputs

def ceil(x: ndarray, /, xp, **kwargs) -> ndarray:
    if xp.issubdtype(x.dtype, xp.integer):
        return x
    return xp.ceil(x, **kwargs)

def floor(x: ndarray, /, xp, **kwargs) -> ndarray:
    if xp.issubdtype(x.dtype, xp.integer):
        return x
    return xp.floor(x, **kwargs)

def trunc(x: ndarray, /, xp, **kwargs) -> ndarray:
    if xp.issubdtype(x.dtype, xp.integer):
        return x
    return xp.trunc(x, **kwargs)

# linear algebra functions

def matmul(x1: ndarray, x2: ndarray, /, xp, **kwargs) -> ndarray:
    return xp.matmul(x1, x2, **kwargs)

# Unlike transpose, matrix_transpose only transposes the last two axes.
def matrix_transpose(x: ndarray, /, xp) -> ndarray:
    if x.ndim < 2:
        raise ValueError("x must be at least 2-dimensional for matrix_transpose")
    return xp.swapaxes(x, -1, -2)

def tensordot(x1: ndarray,
              x2: ndarray,
              /,
              xp,
              *,
              axes: Union[int, Tuple[Sequence[int], Sequence[int]]] = 2,
              **kwargs,
) -> ndarray:
    return xp.tensordot(x1, x2, axes=axes, **kwargs)

def vecdot(x1: ndarray, x2: ndarray, /, xp, *, axis: int = -1) -> ndarray:
    if x1.shape[axis] != x2.shape[axis]:
        raise ValueError("x1 and x2 must have the same size along the given axis")

    if hasattr(xp, 'broadcast_tensors'):
        _broadcast = xp.broadcast_tensors
    else:
        _broadcast = xp.broadcast_arrays

    x1_ = xp.moveaxis(x1, axis, -1)
    x2_ = xp.moveaxis(x2, axis, -1)
    x1_, x2_ = _broadcast(x1_, x2_)

    res = x1_[..., None, :] @ x2_[..., None]
    return res[..., 0, 0]

# isdtype is a new function in the 2022.12 array API specification.

def isdtype(
    dtype: Dtype, kind: Union[Dtype, str, Tuple[Union[Dtype, str], ...]], xp,
    *, _tuple=True, # Disallow nested tuples
) -> bool:
    """
    Returns a boolean indicating whether a provided dtype is of a specified data type ``kind``.

    Note that outside of this function, this compat library does not yet fully
    support complex numbers.

    See
    https://data-apis.org/array-api/latest/API_specification/generated/array_api.isdtype.html
    for more details
    """
    if isinstance(kind, tuple) and _tuple:
        return any(isdtype(dtype, k, xp, _tuple=False) for k in kind)
    elif isinstance(kind, str):
        if kind == 'bool':
            return dtype == xp.bool_
        elif kind == 'signed integer':
            return xp.issubdtype(dtype, xp.signedinteger)
        elif kind == 'unsigned integer':
            return xp.issubdtype(dtype, xp.unsignedinteger)
        elif kind == 'integral':
            return xp.issubdtype(dtype, xp.integer)
        elif kind == 'real floating':
            return xp.issubdtype(dtype, xp.floating)
        elif kind == 'complex floating':
            return xp.issubdtype(dtype, xp.complexfloating)
        elif kind == 'numeric':
            return xp.issubdtype(dtype, xp.number)
        else:
            raise ValueError(f"Unrecognized data type kind: {kind!r}")
    else:
        # This will allow things that aren't required by the spec, like
        # isdtype(np.float64, float) or isdtype(np.int64, 'l'). Should we be
        # more strict here to match the type annotation? Note that the
        # array_api_strict implementation will be very strict.
        return dtype == kind

# unstack is a new function in the 2023.12 array API standard
def unstack(x: ndarray, /, xp, *, axis: int = 0) -> Tuple[ndarray, ...]:
    if x.ndim == 0:
        raise ValueError("Input array must be at least 1-d.")
    return tuple(xp.moveaxis(x, axis, 0))

# numpy 1.26 does not use the standard definition for sign on complex numbers

def sign(x: ndarray, /, xp, **kwargs) -> ndarray:
    if isdtype(x.dtype, 'complex floating', xp=xp):
        out = (x/xp.abs(x, **kwargs))[...]
        # sign(0) = 0 but the above formula would give nan
        out[x == 0+0j] = 0+0j
    else:
        out = xp.sign(x, **kwargs)
    # CuPy sign() does not propagate nans. See
    # https://github.com/data-apis/array-api-compat/issues/136
    if is_cupy_namespace(xp) and isdtype(x.dtype, 'real floating', xp=xp):
        out[xp.isnan(x)] = xp.nan
    return out[()]

__all__ = ['arange', 'empty', 'empty_like', 'eye', 'full', 'full_like',
           'linspace', 'ones', 'ones_like', 'zeros', 'zeros_like',
           'UniqueAllResult', 'UniqueCountsResult', 'UniqueInverseResult',
           'unique_all', 'unique_counts', 'unique_inverse', 'unique_values',
           'astype', 'std', 'var', 'cumulative_sum', 'clip', 'permute_dims',
           'reshape', 'argsort', 'sort', 'nonzero', 'ceil', 'floor', 'trunc',
           'matmul', 'matrix_transpose', 'tensordot', 'vecdot', 'isdtype',
           'unstack', 'sign']
