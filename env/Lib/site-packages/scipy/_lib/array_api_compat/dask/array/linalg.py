from __future__ import annotations

from ...common import _linalg
from ..._internal import get_xp

# Exports
from dask.array.linalg import * # noqa: F403
from dask.array import outer

# These functions are in both the main and linalg namespaces
from dask.array import matmul, tensordot
from ._aliases import matrix_transpose, vecdot

import dask.array as da

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ...common._typing import Array
    from typing import Literal

# dask.array.linalg doesn't have __all__. If it is added, replace this with
#
# from dask.array.linalg import __all__ as linalg_all
_n = {}
exec('from dask.array.linalg import *', _n)
del _n['__builtins__']
if 'annotations' in _n:
    del _n['annotations']
linalg_all = list(_n)
del _n

EighResult = _linalg.EighResult
QRResult = _linalg.QRResult
SlogdetResult = _linalg.SlogdetResult
SVDResult = _linalg.SVDResult
# TODO: use the QR wrapper once dask
# supports the mode keyword on QR
# https://github.com/dask/dask/issues/10388
#qr = get_xp(da)(_linalg.qr)
def qr(x: Array, mode: Literal['reduced', 'complete'] = 'reduced',
       **kwargs) -> QRResult:
    if mode != "reduced":
        raise ValueError("dask arrays only support using mode='reduced'")
    return QRResult(*da.linalg.qr(x, **kwargs))
trace = get_xp(da)(_linalg.trace)
cholesky = get_xp(da)(_linalg.cholesky)
matrix_rank = get_xp(da)(_linalg.matrix_rank)
matrix_norm = get_xp(da)(_linalg.matrix_norm)


# Wrap the svd functions to not pass full_matrices to dask
# when full_matrices=False (as that is the default behavior for dask),
# and dask doesn't have the full_matrices keyword
def svd(x: Array, full_matrices: bool = True, **kwargs) -> SVDResult:
    if full_matrices:
        raise ValueError("full_matrics=True is not supported by dask.")
    return da.linalg.svd(x, coerce_signs=False, **kwargs)

def svdvals(x: Array) -> Array:
    # TODO: can't avoid computing U or V for dask
    _, s, _ =  svd(x)
    return s

vector_norm = get_xp(da)(_linalg.vector_norm)
diagonal = get_xp(da)(_linalg.diagonal)

__all__ = linalg_all + ["trace", "outer", "matmul", "tensordot",
                        "matrix_transpose", "vecdot", "EighResult",
                        "QRResult", "SlogdetResult", "SVDResult", "qr",
                        "cholesky", "matrix_rank", "matrix_norm", "svdvals",
                        "vector_norm", "diagonal"]

_all_ignore = ['get_xp', 'da', 'linalg_all']
