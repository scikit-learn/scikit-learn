import sys
import types
from typing import (
    Literal as L,
    Type,
    Union,
    Tuple,
    overload,
    Any,
    TypeVar,
    Dict,
    List,
    Iterable,
    Protocol,
    TypedDict,
)

from numpy import (
    ndarray,
    dtype,
    generic,
    bool_,
    ubyte,
    ushort,
    uintc,
    uint,
    ulonglong,
    byte,
    short,
    intc,
    int_,
    longlong,
    half,
    single,
    double,
    longdouble,
    csingle,
    cdouble,
    clongdouble,
    datetime64,
    timedelta64,
    object_,
    str_,
    bytes_,
    void,
)

from numpy.core._type_aliases import (
    sctypeDict as sctypeDict,
    sctypes as sctypes,
)

from numpy.typing import DTypeLike, ArrayLike, _SupportsDType

_T = TypeVar("_T")
_SCT = TypeVar("_SCT", bound=generic)

# A paramtrizable subset of `npt.DTypeLike`
_DTypeLike = Union[
    Type[_SCT],
    dtype[_SCT],
    _SupportsDType[dtype[_SCT]],
]

class _CastFunc(Protocol):
    def __call__(
        self, x: ArrayLike, k: DTypeLike = ...
    ) -> ndarray[Any, dtype[Any]]: ...

class _TypeCodes(TypedDict):
    Character: L['c']
    Integer: L['bhilqp']
    UnsignedInteger: L['BHILQP']
    Float: L['efdg']
    Complex: L['FDG']
    AllInteger: L['bBhHiIlLqQpP']
    AllFloat: L['efdgFDG']
    Datetime: L['Mm']
    All: L['?bhilqpBHILQPefdgFDGSUVOMm']

class _typedict(Dict[Type[generic], _T]):
    def __getitem__(self, key: DTypeLike) -> _T: ...

if sys.version_info >= (3, 10):
    _TypeTuple = Union[
        Type[Any],
        types.UnionType,
        Tuple[Union[Type[Any], types.UnionType, Tuple[Any, ...]], ...],
    ]
else:
    _TypeTuple = Union[
        Type[Any],
        Tuple[Union[Type[Any], Tuple[Any, ...]], ...],
    ]

__all__: List[str]

@overload
def maximum_sctype(t: _DTypeLike[_SCT]) -> Type[_SCT]: ...
@overload
def maximum_sctype(t: DTypeLike) -> Type[Any]: ...

@overload
def issctype(rep: dtype[Any] | Type[Any]) -> bool: ...
@overload
def issctype(rep: object) -> L[False]: ...

@overload
def obj2sctype(rep: _DTypeLike[_SCT], default: None = ...) -> None | Type[_SCT]: ...
@overload
def obj2sctype(rep: _DTypeLike[_SCT], default: _T) -> _T | Type[_SCT]: ...
@overload
def obj2sctype(rep: DTypeLike, default: None = ...) -> None | Type[Any]: ...
@overload
def obj2sctype(rep: DTypeLike, default: _T) -> _T | Type[Any]: ...
@overload
def obj2sctype(rep: object, default: None = ...) -> None: ...
@overload
def obj2sctype(rep: object, default: _T) -> _T: ...

@overload
def issubclass_(arg1: Type[Any], arg2: _TypeTuple) -> bool: ...
@overload
def issubclass_(arg1: object, arg2: object) -> L[False]: ...

def issubsctype(arg1: DTypeLike, arg2: DTypeLike) -> bool: ...

def issubdtype(arg1: DTypeLike, arg2: DTypeLike) -> bool: ...

def sctype2char(sctype: DTypeLike) -> str: ...

def find_common_type(
    array_types: Iterable[DTypeLike],
    scalar_types: Iterable[DTypeLike],
) -> dtype[Any]: ...

cast: _typedict[_CastFunc]
nbytes: _typedict[int]
typecodes: _TypeCodes
ScalarType: Tuple[
    Type[int],
    Type[float],
    Type[complex],
    Type[int],
    Type[bool],
    Type[bytes],
    Type[str],
    Type[memoryview],
    Type[bool_],
    Type[csingle],
    Type[cdouble],
    Type[clongdouble],
    Type[half],
    Type[single],
    Type[double],
    Type[longdouble],
    Type[byte],
    Type[short],
    Type[intc],
    Type[int_],
    Type[longlong],
    Type[timedelta64],
    Type[datetime64],
    Type[object_],
    Type[bytes_],
    Type[str_],
    Type[ubyte],
    Type[ushort],
    Type[uintc],
    Type[uint],
    Type[ulonglong],
    Type[void],
]
