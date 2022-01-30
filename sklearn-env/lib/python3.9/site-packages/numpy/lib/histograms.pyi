from typing import (
    Literal as L,
    List,
    Tuple,
    Any,
    SupportsIndex,
    Sequence,
)

from numpy.typing import (
    NDArray,
    ArrayLike,
)

_BinKind = L[
    "stone",
    "auto",
    "doane",
    "fd",
    "rice",
    "scott",
    "sqrt",
    "sturges",
]

__all__: List[str]

def histogram_bin_edges(
    a: ArrayLike,
    bins: _BinKind | SupportsIndex | ArrayLike = ...,
    range: None | Tuple[float, float] = ...,
    weights: None | ArrayLike = ...,
) -> NDArray[Any]: ...

def histogram(
    a: ArrayLike,
    bins: _BinKind | SupportsIndex | ArrayLike = ...,
    range: None | Tuple[float, float] = ...,
    normed: None = ...,
    weights: None | ArrayLike = ...,
    density: bool = ...,
) -> Tuple[NDArray[Any], NDArray[Any]]: ...

def histogramdd(
    sample: ArrayLike,
    bins: SupportsIndex | ArrayLike = ...,
    range: Sequence[Tuple[float, float]] = ...,
    normed: None | bool = ...,
    weights: None | ArrayLike = ...,
    density: None | bool = ...,
) -> Tuple[NDArray[Any], List[NDArray[Any]]]: ...
