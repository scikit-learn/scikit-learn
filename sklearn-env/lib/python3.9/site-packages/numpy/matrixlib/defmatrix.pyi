from typing import List, Any, Sequence, Mapping
from numpy import matrix as matrix
from numpy.typing import ArrayLike, DTypeLike, NDArray

__all__: List[str]

def bmat(
    obj: str | Sequence[ArrayLike] | NDArray[Any],
    ldict: None | Mapping[str, Any] = ...,
    gdict: None | Mapping[str, Any] = ...,
) -> matrix[Any, Any]: ...

def asmatrix(data: ArrayLike, dtype: DTypeLike = ...) -> matrix[Any, Any]: ...

mat = asmatrix
