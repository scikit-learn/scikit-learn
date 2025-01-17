from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple
if TYPE_CHECKING:
    from typing import Literal, Optional, Tuple, Union
    from ._typing import ndarray

import math

import numpy as np
if np.__version__[0] == "2":
    from numpy.lib.array_utils import normalize_axis_tuple
else:
    from numpy.core.numeric import normalize_axis_tuple

from ._aliases import matmul, matrix_transpose, tensordot, vecdot, isdtype
from .._internal import get_xp

# These are in the main NumPy namespace but not in numpy.linalg
def cross(x1: ndarray, x2: ndarray, /, xp, *, axis: int = -1, **kwargs) -> ndarray:
    return xp.cross(x1, x2, axis=axis, **kwargs)

def outer(x1: ndarray, x2: ndarray, /, xp, **kwargs) -> ndarray:
    return xp.outer(x1, x2, **kwargs)

class EighResult(NamedTuple):
    eigenvalues: ndarray
    eigenvectors: ndarray

class QRResult(NamedTuple):
    Q: ndarray
    R: ndarray

class SlogdetResult(NamedTuple):
    sign: ndarray
    logabsdet: ndarray

class SVDResult(NamedTuple):
    U: ndarray
    S: ndarray
    Vh: ndarray

# These functions are the same as their NumPy counterparts except they return
# a namedtuple.
def eigh(x: ndarray, /, xp, **kwargs) -> EighResult:
    return EighResult(*xp.linalg.eigh(x, **kwargs))

def qr(x: ndarray, /, xp, *, mode: Literal['reduced', 'complete'] = 'reduced',
       **kwargs) -> QRResult:
    return QRResult(*xp.linalg.qr(x, mode=mode, **kwargs))

def slogdet(x: ndarray, /, xp, **kwargs) -> SlogdetResult:
    return SlogdetResult(*xp.linalg.slogdet(x, **kwargs))

def svd(x: ndarray, /, xp, *, full_matrices: bool = True, **kwargs) -> SVDResult:
    return SVDResult(*xp.linalg.svd(x, full_matrices=full_matrices, **kwargs))

# These functions have additional keyword arguments

# The upper keyword argument is new from NumPy
def cholesky(x: ndarray, /, xp, *, upper: bool = False, **kwargs) -> ndarray:
    L = xp.linalg.cholesky(x, **kwargs)
    if upper:
        U = get_xp(xp)(matrix_transpose)(L)
        if get_xp(xp)(isdtype)(U.dtype, 'complex floating'):
            U = xp.conj(U)
        return U
    return L

# The rtol keyword argument of matrix_rank() and pinv() is new from NumPy.
# Note that it has a different semantic meaning from tol and rcond.
def matrix_rank(x: ndarray,
                /,
                xp,
                *,
                rtol: Optional[Union[float, ndarray]] = None,
                **kwargs) -> ndarray:
    # this is different from xp.linalg.matrix_rank, which supports 1
    # dimensional arrays.
    if x.ndim < 2:
        raise xp.linalg.LinAlgError("1-dimensional array given. Array must be at least two-dimensional")
    S = get_xp(xp)(svdvals)(x, **kwargs)
    if rtol is None:
        tol = S.max(axis=-1, keepdims=True) * max(x.shape[-2:]) * xp.finfo(S.dtype).eps
    else:
        # this is different from xp.linalg.matrix_rank, which does not
        # multiply the tolerance by the largest singular value.
        tol = S.max(axis=-1, keepdims=True)*xp.asarray(rtol)[..., xp.newaxis]
    return xp.count_nonzero(S > tol, axis=-1)

def pinv(x: ndarray, /, xp, *, rtol: Optional[Union[float, ndarray]] = None, **kwargs) -> ndarray:
    # this is different from xp.linalg.pinv, which does not multiply the
    # default tolerance by max(M, N).
    if rtol is None:
        rtol = max(x.shape[-2:]) * xp.finfo(x.dtype).eps
    return xp.linalg.pinv(x, rcond=rtol, **kwargs)

# These functions are new in the array API spec

def matrix_norm(x: ndarray, /, xp, *, keepdims: bool = False, ord: Optional[Union[int, float, Literal['fro', 'nuc']]] = 'fro') -> ndarray:
    return xp.linalg.norm(x, axis=(-2, -1), keepdims=keepdims, ord=ord)

# svdvals is not in NumPy (but it is in SciPy). It is equivalent to
# xp.linalg.svd(compute_uv=False).
def svdvals(x: ndarray, /, xp) -> Union[ndarray, Tuple[ndarray, ...]]:
    return xp.linalg.svd(x, compute_uv=False)

def vector_norm(x: ndarray, /, xp, *, axis: Optional[Union[int, Tuple[int, ...]]] = None, keepdims: bool = False, ord: Optional[Union[int, float]] = 2) -> ndarray:
    # xp.linalg.norm tries to do a matrix norm whenever axis is a 2-tuple or
    # when axis=None and the input is 2-D, so to force a vector norm, we make
    # it so the input is 1-D (for axis=None), or reshape so that norm is done
    # on a single dimension.
    if axis is None:
        # Note: xp.linalg.norm() doesn't handle 0-D arrays
        _x = x.ravel()
        _axis = 0
    elif isinstance(axis, tuple):
        # Note: The axis argument supports any number of axes, whereas
        # xp.linalg.norm() only supports a single axis for vector norm.
        normalized_axis = normalize_axis_tuple(axis, x.ndim)
        rest = tuple(i for i in range(x.ndim) if i not in normalized_axis)
        newshape = axis + rest
        _x = xp.transpose(x, newshape).reshape(
            (math.prod([x.shape[i] for i in axis]), *[x.shape[i] for i in rest]))
        _axis = 0
    else:
        _x = x
        _axis = axis

    res = xp.linalg.norm(_x, axis=_axis, ord=ord)

    if keepdims:
        # We can't reuse xp.linalg.norm(keepdims) because of the reshape hacks
        # above to avoid matrix norm logic.
        shape = list(x.shape)
        _axis = normalize_axis_tuple(range(x.ndim) if axis is None else axis, x.ndim)
        for i in _axis:
            shape[i] = 1
        res = xp.reshape(res, tuple(shape))

    return res

# xp.diagonal and xp.trace operate on the first two axes whereas these
# operates on the last two

def diagonal(x: ndarray, /, xp, *, offset: int = 0, **kwargs) -> ndarray:
    return xp.diagonal(x, offset=offset, axis1=-2, axis2=-1, **kwargs)

def trace(x: ndarray, /, xp, *, offset: int = 0, dtype=None, **kwargs) -> ndarray:
    return xp.asarray(xp.trace(x, offset=offset, dtype=dtype, axis1=-2, axis2=-1, **kwargs))

__all__ = ['cross', 'matmul', 'outer', 'tensordot', 'EighResult',
           'QRResult', 'SlogdetResult', 'SVDResult', 'eigh', 'qr', 'slogdet',
           'svd', 'cholesky', 'matrix_rank', 'pinv', 'matrix_norm',
           'matrix_transpose', 'svdvals', 'vecdot', 'vector_norm', 'diagonal',
           'trace']
