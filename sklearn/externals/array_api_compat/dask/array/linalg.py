from __future__ import annotations

from typing import Literal

import dask.array as da

# The `matmul` and `tensordot` functions are in both the main and linalg namespaces
from dask.array import matmul, outer, tensordot

# Exports
from ..._internal import clone_module, get_xp
from ...common import _linalg
from ...common._typing import Array

__all__ = clone_module("dask.array.linalg", globals())

from ._aliases import matrix_transpose, vecdot

EighResult = _linalg.EighResult
QRResult = _linalg.QRResult
SlogdetResult = _linalg.SlogdetResult
SVDResult = _linalg.SVDResult
# TODO: use the QR wrapper once dask
# supports the mode keyword on QR
# https://github.com/dask/dask/issues/10388
#qr = get_xp(da)(_linalg.qr)
def qr(  # type: ignore[no-redef]
    x: Array,
    mode: Literal["reduced", "complete"] = "reduced",
    **kwargs: object,
) -> QRResult:
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
def svd(x: Array, full_matrices: bool = True, **kwargs: object) -> SVDResult:  # type: ignore[no-redef]
    if full_matrices:
        raise ValueError("full_matrics=True is not supported by dask.")
    return da.linalg.svd(x, coerce_signs=False, **kwargs)

def svdvals(x: Array) -> Array:
    # TODO: can't avoid computing U or V for dask
    _, s, _ =  svd(x)
    return s

vector_norm = get_xp(da)(_linalg.vector_norm)
diagonal = get_xp(da)(_linalg.diagonal)

__all__ += ["trace", "outer", "matmul", "tensordot",
            "matrix_transpose", "vecdot", "EighResult",
            "QRResult", "SlogdetResult", "SVDResult", "qr",
            "cholesky", "matrix_rank", "matrix_norm", "svdvals",
            "vector_norm", "diagonal"]

def __dir__() -> list[str]:
    return __all__
