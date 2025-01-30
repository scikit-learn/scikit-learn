from collections.abc import Sequence, Mapping
from typing import Any

from numpy import matrix
from numpy._typing import ArrayLike, DTypeLike, NDArray

__all__ = ["matrix", "bmat", "asmatrix"]

def bmat(
    obj: str | Sequence[ArrayLike] | NDArray[Any],
    ldict: None | Mapping[str, Any] = ...,
    gdict: None | Mapping[str, Any] = ...,
) -> matrix[tuple[int, int], Any]: ...

def asmatrix(data: ArrayLike, dtype: DTypeLike = ...) -> matrix[tuple[int, int], Any]: ...

mat = asmatrix
