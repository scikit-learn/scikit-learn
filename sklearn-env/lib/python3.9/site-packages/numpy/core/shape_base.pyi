from typing import TypeVar, overload, List, Sequence, Any, SupportsIndex

from numpy import generic, dtype
from numpy.typing import ArrayLike, NDArray, _FiniteNestedSequence, _SupportsArray

_SCT = TypeVar("_SCT", bound=generic)
_ArrayType = TypeVar("_ArrayType", bound=NDArray[Any])

_ArrayLike = _FiniteNestedSequence[_SupportsArray[dtype[_SCT]]]

__all__: List[str]

@overload
def atleast_1d(arys: _ArrayLike[_SCT], /) -> NDArray[_SCT]: ...
@overload
def atleast_1d(arys: ArrayLike, /) -> NDArray[Any]: ...
@overload
def atleast_1d(*arys: ArrayLike) -> List[NDArray[Any]]: ...

@overload
def atleast_2d(arys: _ArrayLike[_SCT], /) -> NDArray[_SCT]: ...
@overload
def atleast_2d(arys: ArrayLike, /) -> NDArray[Any]: ...
@overload
def atleast_2d(*arys: ArrayLike) -> List[NDArray[Any]]: ...

@overload
def atleast_3d(arys: _ArrayLike[_SCT], /) -> NDArray[_SCT]: ...
@overload
def atleast_3d(arys: ArrayLike, /) -> NDArray[Any]: ...
@overload
def atleast_3d(*arys: ArrayLike) -> List[NDArray[Any]]: ...

@overload
def vstack(tup: Sequence[_ArrayLike[_SCT]]) -> NDArray[_SCT]: ...
@overload
def vstack(tup: Sequence[ArrayLike]) -> NDArray[Any]: ...

@overload
def hstack(tup: Sequence[_ArrayLike[_SCT]]) -> NDArray[_SCT]: ...
@overload
def hstack(tup: Sequence[ArrayLike]) -> NDArray[Any]: ...

@overload
def stack(
    arrays: Sequence[_ArrayLike[_SCT]],
    axis: SupportsIndex = ...,
    out: None = ...,
) -> NDArray[_SCT]: ...
@overload
def stack(
    arrays: Sequence[ArrayLike],
    axis: SupportsIndex = ...,
    out: None = ...,
) -> NDArray[Any]: ...
@overload
def stack(
    arrays: Sequence[ArrayLike],
    axis: SupportsIndex = ...,
    out: _ArrayType = ...,
) -> _ArrayType: ...

@overload
def block(arrays: _ArrayLike[_SCT]) -> NDArray[_SCT]: ...
@overload
def block(arrays: ArrayLike) -> NDArray[Any]: ...
