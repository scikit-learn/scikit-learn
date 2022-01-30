from typing import List, Any, TypeVar, overload

from numpy import generic, dtype, integer, floating, complexfloating
from numpy.typing import (
    NDArray,
    ArrayLike,
    _ShapeLike,
    _SupportsArray,
    _FiniteNestedSequence,
    _ArrayLikeFloat_co,
    _ArrayLikeComplex_co,
)

_SCT = TypeVar("_SCT", bound=generic)

_ArrayLike = _FiniteNestedSequence[_SupportsArray[dtype[_SCT]]]

__all__: List[str]

@overload
def fftshift(x: _ArrayLike[_SCT], axes: None | _ShapeLike = ...) -> NDArray[_SCT]: ...
@overload
def fftshift(x: ArrayLike, axes: None | _ShapeLike = ...) -> NDArray[Any]: ...

@overload
def ifftshift(x: _ArrayLike[_SCT], axes: None | _ShapeLike = ...) -> NDArray[_SCT]: ...
@overload
def ifftshift(x: ArrayLike, axes: None | _ShapeLike = ...) -> NDArray[Any]: ...

@overload
def fftfreq(
    n: int | integer[Any],
    d: _ArrayLikeFloat_co,
) -> NDArray[floating[Any]]: ...
@overload
def fftfreq(
    n: int | integer[Any],
    d: _ArrayLikeComplex_co,
) -> NDArray[complexfloating[Any, Any]]: ...

@overload
def rfftfreq(
    n: int | integer[Any],
    d: _ArrayLikeFloat_co,
) -> NDArray[floating[Any]]: ...
@overload
def rfftfreq(
    n: int | integer[Any],
    d: _ArrayLikeComplex_co,
) -> NDArray[complexfloating[Any, Any]]: ...
