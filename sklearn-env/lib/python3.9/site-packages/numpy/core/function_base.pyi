from typing import overload, Tuple, Union, Sequence, Any, SupportsIndex, Literal, List

from numpy import ndarray
from numpy.typing import ArrayLike, DTypeLike, _SupportsArray, _NumberLike_co

# TODO: wait for support for recursive types
_ArrayLikeNested = Sequence[Sequence[Any]]
_ArrayLikeNumber = Union[
    _NumberLike_co, Sequence[_NumberLike_co], ndarray, _SupportsArray, _ArrayLikeNested
]

__all__: List[str]

@overload
def linspace(
    start: _ArrayLikeNumber,
    stop: _ArrayLikeNumber,
    num: SupportsIndex = ...,
    endpoint: bool = ...,
    retstep: Literal[False] = ...,
    dtype: DTypeLike = ...,
    axis: SupportsIndex = ...,
) -> ndarray: ...
@overload
def linspace(
    start: _ArrayLikeNumber,
    stop: _ArrayLikeNumber,
    num: SupportsIndex = ...,
    endpoint: bool = ...,
    retstep: Literal[True] = ...,
    dtype: DTypeLike = ...,
    axis: SupportsIndex = ...,
) -> Tuple[ndarray, Any]: ...

def logspace(
    start: _ArrayLikeNumber,
    stop: _ArrayLikeNumber,
    num: SupportsIndex = ...,
    endpoint: bool = ...,
    base: _ArrayLikeNumber = ...,
    dtype: DTypeLike = ...,
    axis: SupportsIndex = ...,
) -> ndarray: ...

def geomspace(
    start: _ArrayLikeNumber,
    stop: _ArrayLikeNumber,
    num: SupportsIndex = ...,
    endpoint: bool = ...,
    dtype: DTypeLike = ...,
    axis: SupportsIndex = ...,
) -> ndarray: ...

# Re-exported to `np.lib.function_base`
def add_newdoc(
    place: str,
    obj: str,
    doc: str | Tuple[str, str] | List[Tuple[str, str]],
    warn_on_python: bool = ...,
) -> None: ...
