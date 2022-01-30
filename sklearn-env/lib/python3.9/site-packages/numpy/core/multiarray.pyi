# TODO: Sort out any and all missing functions in this namespace

import os
import datetime as dt
from typing import (
    Literal as L,
    Any,
    Callable,
    Iterable,
    Optional,
    overload,
    TypeVar,
    List,
    Type,
    Union,
    Sequence,
    Tuple,
    SupportsIndex,
    final,
    Final,
    Protocol,
)

from numpy import (
    # Re-exports
    busdaycalendar as busdaycalendar,
    broadcast as broadcast,
    dtype as dtype,
    ndarray as ndarray,
    nditer as nditer,

    # The rest
    ufunc,
    str_,
    bool_,
    uint8,
    intp,
    int_,
    float64,
    timedelta64,
    datetime64,
    generic,
    unsignedinteger,
    signedinteger,
    floating,
    complexfloating,
    _OrderKACF,
    _OrderCF,
    _CastingKind,
    _ModeKind,
    _SupportsBuffer,
    _IOProtocol,
    _CopyMode,
    _NDIterFlagsKind,
    _NDIterOpFlagsKind,
)

from numpy.typing import (
    # Shapes
    _ShapeLike,

    # DTypes
    DTypeLike,
    _SupportsDType,

    # Arrays
    NDArray,
    ArrayLike,
    _SupportsArray,
    _NestedSequence,
    _FiniteNestedSequence,
    _ArrayLikeBool_co,
    _ArrayLikeUInt_co,
    _ArrayLikeInt_co,
    _ArrayLikeFloat_co,
    _ArrayLikeComplex_co,
    _ArrayLikeTD64_co,
    _ArrayLikeDT64_co,
    _ArrayLikeObject_co,
    _ArrayLikeStr_co,
    _ArrayLikeBytes_co,
    _ScalarLike_co,
    _IntLike_co,
    _FloatLike_co,
    _TD64Like_co,
)

_SCT = TypeVar("_SCT", bound=generic)
_ArrayType = TypeVar("_ArrayType", bound=NDArray[Any])

# Subscriptable subsets of `npt.DTypeLike` and `npt.ArrayLike`
_DTypeLike = Union[
    dtype[_SCT],
    Type[_SCT],
    _SupportsDType[dtype[_SCT]],
]
_ArrayLike = _FiniteNestedSequence[_SupportsArray[dtype[_SCT]]]

# Valid time units
_UnitKind = L[
    "Y",
    "M",
    "D",
    "h",
    "m",
    "s",
    "ms",
    "us", "μs",
    "ns",
    "ps",
    "fs",
    "as",
]
_RollKind = L[  # `raise` is deliberately excluded
    "nat",
    "forward",
    "following",
    "backward",
    "preceding",
    "modifiedfollowing",
    "modifiedpreceding",
]

__all__: List[str]

ALLOW_THREADS: Final[int]  # 0 or 1 (system-specific)
BUFSIZE: L[8192]
CLIP: L[0]
WRAP: L[1]
RAISE: L[2]
MAXDIMS: L[32]
MAY_SHARE_BOUNDS: L[0]
MAY_SHARE_EXACT: L[-1]
tracemalloc_domain: L[389047]

@overload
def empty_like(
    prototype: _ArrayType,
    dtype: None = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    shape: Optional[_ShapeLike] = ...,
) -> _ArrayType: ...
@overload
def empty_like(
    prototype: _ArrayLike[_SCT],
    dtype: None = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    shape: Optional[_ShapeLike] = ...,
) -> NDArray[_SCT]: ...
@overload
def empty_like(
    prototype: object,
    dtype: None = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    shape: Optional[_ShapeLike] = ...,
) -> NDArray[Any]: ...
@overload
def empty_like(
    prototype: Any,
    dtype: _DTypeLike[_SCT],
    order: _OrderKACF = ...,
    subok: bool = ...,
    shape: Optional[_ShapeLike] = ...,
) -> NDArray[_SCT]: ...
@overload
def empty_like(
    prototype: Any,
    dtype: DTypeLike,
    order: _OrderKACF = ...,
    subok: bool = ...,
    shape: Optional[_ShapeLike] = ...,
) -> NDArray[Any]: ...

@overload
def array(
    object: _ArrayType,
    dtype: None = ...,
    *,
    copy: bool | _CopyMode = ...,
    order: _OrderKACF = ...,
    subok: L[True],
    ndmin: int = ...,
    like: ArrayLike = ...,
) -> _ArrayType: ...
@overload
def array(
    object: _ArrayLike[_SCT],
    dtype: None = ...,
    *,
    copy: bool | _CopyMode = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    ndmin: int = ...,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def array(
    object: object,
    dtype: None = ...,
    *,
    copy: bool | _CopyMode = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    ndmin: int = ...,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def array(
    object: Any,
    dtype: _DTypeLike[_SCT],
    *,
    copy: bool | _CopyMode = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    ndmin: int = ...,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def array(
    object: Any,
    dtype: DTypeLike,
    *,
    copy: bool | _CopyMode = ...,
    order: _OrderKACF = ...,
    subok: bool = ...,
    ndmin: int = ...,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def zeros(
    shape: _ShapeLike,
    dtype: None = ...,
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[float64]: ...
@overload
def zeros(
    shape: _ShapeLike,
    dtype: _DTypeLike[_SCT],
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def zeros(
    shape: _ShapeLike,
    dtype: DTypeLike,
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def empty(
    shape: _ShapeLike,
    dtype: None = ...,
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[float64]: ...
@overload
def empty(
    shape: _ShapeLike,
    dtype: _DTypeLike[_SCT],
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def empty(
    shape: _ShapeLike,
    dtype: DTypeLike,
    order: _OrderCF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def unravel_index(  # type: ignore[misc]
    indices: _IntLike_co,
    shape: _ShapeLike,
    order: _OrderCF = ...,
) -> Tuple[intp, ...]: ...
@overload
def unravel_index(
    indices: _ArrayLikeInt_co,
    shape: _ShapeLike,
    order: _OrderCF = ...,
) -> Tuple[NDArray[intp], ...]: ...

@overload
def ravel_multi_index(  # type: ignore[misc]
    multi_index: Sequence[_IntLike_co],
    dims: Sequence[SupportsIndex],
    mode: Union[_ModeKind, Tuple[_ModeKind, ...]] = ...,
    order: _OrderCF = ...,
) -> intp: ...
@overload
def ravel_multi_index(
    multi_index: Sequence[_ArrayLikeInt_co],
    dims: Sequence[SupportsIndex],
    mode: Union[_ModeKind, Tuple[_ModeKind, ...]] = ...,
    order: _OrderCF = ...,
) -> NDArray[intp]: ...

@overload
def concatenate(  # type: ignore[misc]
    arrays: _ArrayLike[_SCT],
    /,
    axis: Optional[SupportsIndex] = ...,
    out: None = ...,
    *,
    dtype: None = ...,
    casting: Optional[_CastingKind] = ...
) -> NDArray[_SCT]: ...
@overload
def concatenate(  # type: ignore[misc]
    arrays: ArrayLike,
    /,
    axis: Optional[SupportsIndex] = ...,
    out: None = ...,
    *,
    dtype: None = ...,
    casting: Optional[_CastingKind] = ...
) -> NDArray[Any]: ...
@overload
def concatenate(  # type: ignore[misc]
    arrays: ArrayLike,
    /,
    axis: Optional[SupportsIndex] = ...,
    out: None = ...,
    *,
    dtype: _DTypeLike[_SCT],
    casting: Optional[_CastingKind] = ...
) -> NDArray[_SCT]: ...
@overload
def concatenate(  # type: ignore[misc]
    arrays: ArrayLike,
    /,
    axis: Optional[SupportsIndex] = ...,
    out: None = ...,
    *,
    dtype: DTypeLike,
    casting: Optional[_CastingKind] = ...
) -> NDArray[Any]: ...
@overload
def concatenate(
    arrays: ArrayLike,
    /,
    axis: Optional[SupportsIndex] = ...,
    out: _ArrayType = ...,
    *,
    dtype: DTypeLike = ...,
    casting: Optional[_CastingKind] = ...
) -> _ArrayType: ...

def inner(
    a: ArrayLike,
    b: ArrayLike,
    /,
) -> Any: ...

@overload
def where(
    condition: ArrayLike,
    /,
) -> Tuple[NDArray[intp], ...]: ...
@overload
def where(
    condition: ArrayLike,
    x: ArrayLike,
    y: ArrayLike,
    /,
) -> NDArray[Any]: ...

def lexsort(
    keys: ArrayLike,
    axis: Optional[SupportsIndex] = ...,
) -> Any: ...

def can_cast(
    from_: Union[ArrayLike, DTypeLike],
    to: DTypeLike,
    casting: Optional[_CastingKind] = ...,
) -> bool: ...

def min_scalar_type(
    a: ArrayLike, /,
) -> dtype[Any]: ...

def result_type(
    *arrays_and_dtypes: Union[ArrayLike, DTypeLike],
) -> dtype[Any]: ...

@overload
def dot(a: ArrayLike, b: ArrayLike, out: None = ...) -> Any: ...
@overload
def dot(a: ArrayLike, b: ArrayLike, out: _ArrayType) -> _ArrayType: ...

@overload
def vdot(a: _ArrayLikeBool_co, b: _ArrayLikeBool_co, /) -> bool_: ...  # type: ignore[misc]
@overload
def vdot(a: _ArrayLikeUInt_co, b: _ArrayLikeUInt_co, /) -> unsignedinteger[Any]: ...  # type: ignore[misc]
@overload
def vdot(a: _ArrayLikeInt_co, b: _ArrayLikeInt_co, /) -> signedinteger[Any]: ... # type: ignore[misc]
@overload
def vdot(a: _ArrayLikeFloat_co, b: _ArrayLikeFloat_co, /) -> floating[Any]: ...  # type: ignore[misc]
@overload
def vdot(a: _ArrayLikeComplex_co, b: _ArrayLikeComplex_co, /) -> complexfloating[Any, Any]: ...  # type: ignore[misc]
@overload
def vdot(a: _ArrayLikeTD64_co, b: _ArrayLikeTD64_co, /) -> timedelta64: ...
@overload
def vdot(a: _ArrayLikeObject_co, b: Any, /) -> Any: ...
@overload
def vdot(a: Any, b: _ArrayLikeObject_co, /) -> Any: ...

def bincount(
    x: ArrayLike,
    /,
    weights: Optional[ArrayLike] = ...,
    minlength: SupportsIndex = ...,
) -> NDArray[intp]: ...

def copyto(
    dst: NDArray[Any],
    src: ArrayLike,
    casting: Optional[_CastingKind] = ...,
    where: Optional[_ArrayLikeBool_co] = ...,
) -> None: ...

def putmask(
    a: NDArray[Any],
    mask: _ArrayLikeBool_co,
    values: ArrayLike,
) -> None: ...

def packbits(
    a: _ArrayLikeInt_co,
    /,
    axis: Optional[SupportsIndex] = ...,
    bitorder: L["big", "little"] = ...,
) -> NDArray[uint8]: ...

def unpackbits(
    a: _ArrayLike[uint8],
    /,
    axis: Optional[SupportsIndex] = ...,
    count: Optional[SupportsIndex] = ...,
    bitorder: L["big", "little"] = ...,
) -> NDArray[uint8]: ...

def shares_memory(
    a: object,
    b: object,
    /,
    max_work: Optional[int] = ...,
) -> bool: ...

def may_share_memory(
    a: object,
    b: object,
    /,
    max_work: Optional[int] = ...,
) -> bool: ...

@overload
def asarray(
    a: _ArrayLike[_SCT],
    dtype: None = ...,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asarray(
    a: object,
    dtype: None = ...,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def asarray(
    a: Any,
    dtype: _DTypeLike[_SCT],
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asarray(
    a: Any,
    dtype: DTypeLike,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def asanyarray(
    a: _ArrayType,  # Preserve subclass-information
    dtype: None = ...,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> _ArrayType: ...
@overload
def asanyarray(
    a: _ArrayLike[_SCT],
    dtype: None = ...,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asanyarray(
    a: object,
    dtype: None = ...,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def asanyarray(
    a: Any,
    dtype: _DTypeLike[_SCT],
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asanyarray(
    a: Any,
    dtype: DTypeLike,
    order: _OrderKACF = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def ascontiguousarray(
    a: _ArrayLike[_SCT],
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def ascontiguousarray(
    a: object,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def ascontiguousarray(
    a: Any,
    dtype: _DTypeLike[_SCT],
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def ascontiguousarray(
    a: Any,
    dtype: DTypeLike,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def asfortranarray(
    a: _ArrayLike[_SCT],
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asfortranarray(
    a: object,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def asfortranarray(
    a: Any,
    dtype: _DTypeLike[_SCT],
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def asfortranarray(
    a: Any,
    dtype: DTypeLike,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

# In practice `List[Any]` is list with an int, int and a valid
# `np.seterrcall()` object
def geterrobj() -> List[Any]: ...
def seterrobj(errobj: List[Any], /) -> None: ...

def promote_types(__type1: DTypeLike, __type2: DTypeLike) -> dtype[Any]: ...

# `sep` is a de facto mandatory argument, as its default value is deprecated
@overload
def fromstring(
    string: str | bytes,
    dtype: None = ...,
    count: SupportsIndex = ...,
    *,
    sep: str,
    like: ArrayLike = ...,
) -> NDArray[float64]: ...
@overload
def fromstring(
    string: str | bytes,
    dtype: _DTypeLike[_SCT],
    count: SupportsIndex = ...,
    *,
    sep: str,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def fromstring(
    string: str | bytes,
    dtype: DTypeLike,
    count: SupportsIndex = ...,
    *,
    sep: str,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

def frompyfunc(
    func: Callable[..., Any], /,
    nin: SupportsIndex,
    nout: SupportsIndex,
    *,
    identity: Any = ...,
) -> ufunc: ...

@overload
def fromfile(
    file: str | bytes | os.PathLike[Any] | _IOProtocol,
    dtype: None = ...,
    count: SupportsIndex = ...,
    sep: str = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[float64]: ...
@overload
def fromfile(
    file: str | bytes | os.PathLike[Any] | _IOProtocol,
    dtype: _DTypeLike[_SCT],
    count: SupportsIndex = ...,
    sep: str = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def fromfile(
    file: str | bytes | os.PathLike[Any] | _IOProtocol,
    dtype: DTypeLike,
    count: SupportsIndex = ...,
    sep: str = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def fromiter(
    iter: Iterable[Any],
    dtype: _DTypeLike[_SCT],
    count: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def fromiter(
    iter: Iterable[Any],
    dtype: DTypeLike,
    count: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def frombuffer(
    buffer: _SupportsBuffer,
    dtype: None = ...,
    count: SupportsIndex = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[float64]: ...
@overload
def frombuffer(
    buffer: _SupportsBuffer,
    dtype: _DTypeLike[_SCT],
    count: SupportsIndex = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def frombuffer(
    buffer: _SupportsBuffer,
    dtype: DTypeLike,
    count: SupportsIndex = ...,
    offset: SupportsIndex = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

@overload
def arange(  # type: ignore[misc]
    stop: _IntLike_co,
    /, *,
    dtype: None = ...,
    like: ArrayLike = ...,
) -> NDArray[signedinteger[Any]]: ...
@overload
def arange(  # type: ignore[misc]
    start: _IntLike_co,
    stop: _IntLike_co,
    step: _IntLike_co = ...,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[signedinteger[Any]]: ...
@overload
def arange(  # type: ignore[misc]
    stop: _FloatLike_co,
    /, *,
    dtype: None = ...,
    like: ArrayLike = ...,
) -> NDArray[floating[Any]]: ...
@overload
def arange(  # type: ignore[misc]
    start: _FloatLike_co,
    stop: _FloatLike_co,
    step: _FloatLike_co = ...,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[floating[Any]]: ...
@overload
def arange(
    stop: _TD64Like_co,
    /, *,
    dtype: None = ...,
    like: ArrayLike = ...,
) -> NDArray[timedelta64]: ...
@overload
def arange(
    start: _TD64Like_co,
    stop: _TD64Like_co,
    step: _TD64Like_co = ...,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[timedelta64]: ...
@overload
def arange(  # both start and stop must always be specified for datetime64
    start: datetime64,
    stop: datetime64,
    step: datetime64 = ...,
    dtype: None = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[datetime64]: ...
@overload
def arange(
    stop: Any,
    /, *,
    dtype: _DTypeLike[_SCT],
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def arange(
    start: Any,
    stop: Any,
    step: Any = ...,
    dtype: _DTypeLike[_SCT] = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[_SCT]: ...
@overload
def arange(
    stop: Any, /,
    *,
    dtype: DTypeLike,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...
@overload
def arange(
    start: Any,
    stop: Any,
    step: Any = ...,
    dtype: DTypeLike = ...,
    *,
    like: ArrayLike = ...,
) -> NDArray[Any]: ...

def datetime_data(
    dtype: str | _DTypeLike[datetime64] | _DTypeLike[timedelta64], /,
) -> Tuple[str, int]: ...

# The datetime functions perform unsafe casts to `datetime64[D]`,
# so a lot of different argument types are allowed here

@overload
def busday_count(  # type: ignore[misc]
    begindates: _ScalarLike_co | dt.date,
    enddates: _ScalarLike_co | dt.date,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> int_: ...
@overload
def busday_count(  # type: ignore[misc]
    begindates: ArrayLike | dt.date | _NestedSequence[dt.date],
    enddates: ArrayLike | dt.date | _NestedSequence[dt.date],
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> NDArray[int_]: ...
@overload
def busday_count(
    begindates: ArrayLike | dt.date | _NestedSequence[dt.date],
    enddates: ArrayLike | dt.date | _NestedSequence[dt.date],
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: _ArrayType = ...,
) -> _ArrayType: ...

# `roll="raise"` is (more or less?) equivalent to `casting="safe"`
@overload
def busday_offset(  # type: ignore[misc]
    dates: datetime64 | dt.date,
    offsets: _TD64Like_co | dt.timedelta,
    roll: L["raise"] = ...,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> datetime64: ...
@overload
def busday_offset(  # type: ignore[misc]
    dates: _ArrayLike[datetime64] | dt.date | _NestedSequence[dt.date],
    offsets: _ArrayLikeTD64_co | dt.timedelta | _NestedSequence[dt.timedelta],
    roll: L["raise"] = ...,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> NDArray[datetime64]: ...
@overload
def busday_offset(  # type: ignore[misc]
    dates: _ArrayLike[datetime64] | dt.date | _NestedSequence[dt.date],
    offsets: _ArrayLikeTD64_co | dt.timedelta | _NestedSequence[dt.timedelta],
    roll: L["raise"] = ...,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: _ArrayType = ...,
) -> _ArrayType: ...
@overload
def busday_offset(  # type: ignore[misc]
    dates: _ScalarLike_co | dt.date,
    offsets: _ScalarLike_co | dt.timedelta,
    roll: _RollKind,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> datetime64: ...
@overload
def busday_offset(  # type: ignore[misc]
    dates: ArrayLike | dt.date | _NestedSequence[dt.date],
    offsets: ArrayLike | dt.timedelta | _NestedSequence[dt.timedelta],
    roll: _RollKind,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> NDArray[datetime64]: ...
@overload
def busday_offset(
    dates: ArrayLike | dt.date | _NestedSequence[dt.date],
    offsets: ArrayLike | dt.timedelta | _NestedSequence[dt.timedelta],
    roll: _RollKind,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: _ArrayType = ...,
) -> _ArrayType: ...

@overload
def is_busday(  # type: ignore[misc]
    dates: _ScalarLike_co | dt.date,
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> bool_: ...
@overload
def is_busday(  # type: ignore[misc]
    dates: ArrayLike | _NestedSequence[dt.date],
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: None = ...,
) -> NDArray[bool_]: ...
@overload
def is_busday(
    dates: ArrayLike | _NestedSequence[dt.date],
    weekmask: ArrayLike = ...,
    holidays: None | ArrayLike | dt.date | _NestedSequence[dt.date] = ...,
    busdaycal: None | busdaycalendar = ...,
    out: _ArrayType = ...,
) -> _ArrayType: ...

@overload
def datetime_as_string(  # type: ignore[misc]
    arr: datetime64 | dt.date,
    unit: None | L["auto"] | _UnitKind = ...,
    timezone: L["naive", "UTC", "local"] | dt.tzinfo = ...,
    casting: _CastingKind = ...,
) -> str_: ...
@overload
def datetime_as_string(
    arr: _ArrayLikeDT64_co | _NestedSequence[dt.date],
    unit: None | L["auto"] | _UnitKind = ...,
    timezone: L["naive", "UTC", "local"] | dt.tzinfo = ...,
    casting: _CastingKind = ...,
) -> NDArray[str_]: ...

@overload
def compare_chararrays(
    a1: _ArrayLikeStr_co,
    a2: _ArrayLikeStr_co,
    cmp: L["<", "<=", "==", ">=", ">", "!="],
    rstrip: bool,
) -> NDArray[bool_]: ...
@overload
def compare_chararrays(
    a1: _ArrayLikeBytes_co,
    a2: _ArrayLikeBytes_co,
    cmp: L["<", "<=", "==", ">=", ">", "!="],
    rstrip: bool,
) -> NDArray[bool_]: ...

def add_docstring(obj: Callable[..., Any], docstring: str, /) -> None: ...

_GetItemKeys = L[
    "C", "CONTIGUOUS", "C_CONTIGUOUS",
    "F", "FORTRAN", "F_CONTIGUOUS",
    "W", "WRITEABLE",
    "B", "BEHAVED",
    "O", "OWNDATA",
    "A", "ALIGNED",
    "X", "WRITEBACKIFCOPY",
    "CA", "CARRAY",
    "FA", "FARRAY",
    "FNC",
    "FORC",
]
_SetItemKeys = L[
    "A", "ALIGNED",
    "W", "WRITABLE",
    "X", "WRITEBACKIFCOPY",
]

@final
class flagsobj:
    __hash__: None  # type: ignore[assignment]
    aligned: bool
    # NOTE: deprecated
    # updateifcopy: bool
    writeable: bool
    writebackifcopy: bool
    @property
    def behaved(self) -> bool: ...
    @property
    def c_contiguous(self) -> bool: ...
    @property
    def carray(self) -> bool: ...
    @property
    def contiguous(self) -> bool: ...
    @property
    def f_contiguous(self) -> bool: ...
    @property
    def farray(self) -> bool: ...
    @property
    def fnc(self) -> bool: ...
    @property
    def forc(self) -> bool: ...
    @property
    def fortran(self) -> bool: ...
    @property
    def num(self) -> int: ...
    @property
    def owndata(self) -> bool: ...
    def __getitem__(self, key: _GetItemKeys) -> bool: ...
    def __setitem__(self, key: _SetItemKeys, value: bool) -> None: ...

def nested_iters(
    op: ArrayLike | Sequence[ArrayLike],
    axes: Sequence[Sequence[SupportsIndex]],
    flags: None | Sequence[_NDIterFlagsKind] = ...,
    op_flags: None | Sequence[Sequence[_NDIterOpFlagsKind]] = ...,
    op_dtypes: DTypeLike | Sequence[DTypeLike] = ...,
    order: _OrderKACF = ...,
    casting: _CastingKind = ...,
    buffersize: SupportsIndex = ...,
) -> Tuple[nditer, ...]: ...
