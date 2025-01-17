from collections.abc import Generator
from types import EllipsisType
from typing import (
    Any,
    TypeAlias,
    TypeVar,
    overload,
)

from numpy import ndarray, dtype, generic
from numpy._typing import DTypeLike, NDArray, _Shape as _AnyShape

__all__ = ["Arrayterator"]

# TODO: Rename to ``_ShapeType``
_Shape = TypeVar("_Shape", bound=_AnyShape)
_DType = TypeVar("_DType", bound=dtype[Any])
_ScalarType = TypeVar("_ScalarType", bound=generic)

_Index: TypeAlias = (
    EllipsisType
    | int
    | slice
    | tuple[EllipsisType | int | slice, ...]
)


# NOTE: In reality `Arrayterator` does not actually inherit from `ndarray`,
# but its ``__getattr__` method does wrap around the former and thus has
# access to all its methods

class Arrayterator(ndarray[_Shape, _DType]):
    var: ndarray[_Shape, _DType]  # type: ignore[assignment]
    buf_size: None | int
    start: list[int]
    stop: list[int]
    step: list[int]

    @property  # type: ignore[misc]
    def shape(self) -> tuple[int, ...]: ...
    @property
    def flat(self: NDArray[_ScalarType]) -> Generator[_ScalarType, None, None]: ...
    def __init__(
        self, var: ndarray[_Shape, _DType], buf_size: None | int = ...
    ) -> None: ...
    @overload
    def __array__(self, dtype: None = ..., copy: None | bool = ...) -> ndarray[_AnyShape, _DType]: ...
    @overload
    def __array__(self, dtype: DTypeLike, copy: None | bool = ...) -> NDArray[Any]: ...
    def __getitem__(self, index: _Index) -> Arrayterator[_AnyShape, _DType]: ...
    def __iter__(self) -> Generator[ndarray[_AnyShape, _DType], None, None]: ...
