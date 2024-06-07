from typing import Any, TypeVar, overload

from numpy import bool_, floating, ndarray, object_
from numpy._typing import (
    NDArray,
    _ArrayLikeFloat_co,
    _ArrayLikeObject_co,
    _FloatLike_co,
)

_ArrayType = TypeVar("_ArrayType", bound=ndarray[Any, Any])

__all__: list[str]

@overload
def fix(  # type: ignore[misc]
    x: _FloatLike_co,
    out: None = ...,
) -> floating[Any]: ...
@overload
def fix(
    x: _ArrayLikeFloat_co,
    out: None = ...,
) -> NDArray[floating[Any]]: ...
@overload
def fix(
    x: _ArrayLikeObject_co,
    out: None = ...,
) -> NDArray[object_]: ...
@overload
def fix(
    x: _ArrayLikeFloat_co | _ArrayLikeObject_co,
    out: _ArrayType,
) -> _ArrayType: ...
@overload
def isposinf(  # type: ignore[misc]
    x: _FloatLike_co,
    out: None = ...,
) -> bool_: ...
@overload
def isposinf(
    x: _ArrayLikeFloat_co,
    out: None = ...,
) -> NDArray[bool_]: ...
@overload
def isposinf(
    x: _ArrayLikeFloat_co,
    out: _ArrayType,
) -> _ArrayType: ...
@overload
def isneginf(  # type: ignore[misc]
    x: _FloatLike_co,
    out: None = ...,
) -> bool_: ...
@overload
def isneginf(
    x: _ArrayLikeFloat_co,
    out: None = ...,
) -> NDArray[bool_]: ...
@overload
def isneginf(
    x: _ArrayLikeFloat_co,
    out: _ArrayType,
) -> _ArrayType: ...
