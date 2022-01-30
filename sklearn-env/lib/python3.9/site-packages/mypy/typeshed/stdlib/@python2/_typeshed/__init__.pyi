# Utility types for typeshed

# This module contains various common types to be used by typeshed. The
# module and its types do not exist at runtime. You can use this module
# outside of typeshed, but no API stability guarantees are made. To use
# it in implementation (.py) files, the following construct must be used:
#
#     from typing import TYPE_CHECKING
#     if TYPE_CHECKING:
#         from _typeshed import ...
#
# If on Python versions < 3.10 and "from __future__ import annotations"
# is not used, types from this module must be quoted.

import array
import mmap
from typing import Any, Container, Iterable, Protocol, Text, Tuple, TypeVar, Union
from typing_extensions import Literal, final

_KT = TypeVar("_KT")
_KT_co = TypeVar("_KT_co", covariant=True)
_KT_contra = TypeVar("_KT_contra", contravariant=True)
_VT = TypeVar("_VT")
_VT_co = TypeVar("_VT_co", covariant=True)
_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)
_T_contra = TypeVar("_T_contra", contravariant=True)

# Use for "self" annotations:
#   def __enter__(self: Self) -> Self: ...
Self = TypeVar("Self")  # noqa Y001

class IdentityFunction(Protocol):
    def __call__(self, __x: _T) -> _T: ...

class SupportsLessThan(Protocol):
    def __lt__(self, __other: Any) -> bool: ...

SupportsLessThanT = TypeVar("SupportsLessThanT", bound=SupportsLessThan)  # noqa: Y001

class SupportsDivMod(Protocol[_T_contra, _T_co]):
    def __divmod__(self, __other: _T_contra) -> _T_co: ...

class SupportsRDivMod(Protocol[_T_contra, _T_co]):
    def __rdivmod__(self, __other: _T_contra) -> _T_co: ...

# Mapping-like protocols

class SupportsItems(Protocol[_KT_co, _VT_co]):
    # We want dictionaries to support this on Python 2.
    def items(self) -> Iterable[Tuple[_KT_co, _VT_co]]: ...

class SupportsKeysAndGetItem(Protocol[_KT, _VT_co]):
    def keys(self) -> Iterable[_KT]: ...
    def __getitem__(self, __k: _KT) -> _VT_co: ...

class SupportsGetItem(Container[_KT_contra], Protocol[_KT_contra, _VT_co]):
    def __getitem__(self, __k: _KT_contra) -> _VT_co: ...

class SupportsItemAccess(SupportsGetItem[_KT_contra, _VT], Protocol[_KT_contra, _VT]):
    def __setitem__(self, __k: _KT_contra, __v: _VT) -> None: ...
    def __delitem__(self, __v: _KT_contra) -> None: ...

# These aliases can be used in places where a PathLike object can be used
# instead of a string in Python 3.
StrPath = Text
BytesPath = str
StrOrBytesPath = Text
AnyPath = StrOrBytesPath  # obsolete, will be removed soon

OpenTextModeUpdating = Literal[
    "r+",
    "+r",
    "rt+",
    "r+t",
    "+rt",
    "tr+",
    "t+r",
    "+tr",
    "w+",
    "+w",
    "wt+",
    "w+t",
    "+wt",
    "tw+",
    "t+w",
    "+tw",
    "a+",
    "+a",
    "at+",
    "a+t",
    "+at",
    "ta+",
    "t+a",
    "+ta",
    "x+",
    "+x",
    "xt+",
    "x+t",
    "+xt",
    "tx+",
    "t+x",
    "+tx",
]
OpenTextModeWriting = Literal["w", "wt", "tw", "a", "at", "ta", "x", "xt", "tx"]
OpenTextModeReading = Literal["r", "rt", "tr", "U", "rU", "Ur", "rtU", "rUt", "Urt", "trU", "tUr", "Utr"]
OpenTextMode = Union[OpenTextModeUpdating, OpenTextModeWriting, OpenTextModeReading]
OpenBinaryModeUpdating = Literal[
    "rb+",
    "r+b",
    "+rb",
    "br+",
    "b+r",
    "+br",
    "wb+",
    "w+b",
    "+wb",
    "bw+",
    "b+w",
    "+bw",
    "ab+",
    "a+b",
    "+ab",
    "ba+",
    "b+a",
    "+ba",
    "xb+",
    "x+b",
    "+xb",
    "bx+",
    "b+x",
    "+bx",
]
OpenBinaryModeWriting = Literal["wb", "bw", "ab", "ba", "xb", "bx"]
OpenBinaryModeReading = Literal["rb", "br", "rbU", "rUb", "Urb", "brU", "bUr", "Ubr"]
OpenBinaryMode = Union[OpenBinaryModeUpdating, OpenBinaryModeReading, OpenBinaryModeWriting]

class HasFileno(Protocol):
    def fileno(self) -> int: ...

FileDescriptor = int
FileDescriptorLike = Union[int, HasFileno]

class SupportsRead(Protocol[_T_co]):
    def read(self, __length: int = ...) -> _T_co: ...

class SupportsReadline(Protocol[_T_co]):
    def readline(self, __length: int = ...) -> _T_co: ...

class SupportsNoArgReadline(Protocol[_T_co]):
    def readline(self) -> _T_co: ...

class SupportsWrite(Protocol[_T_contra]):
    def write(self, __s: _T_contra) -> Any: ...

ReadableBuffer = Union[bytes, bytearray, memoryview, array.array[Any], mmap.mmap, buffer]
WriteableBuffer = Union[bytearray, memoryview, array.array[Any], mmap.mmap, buffer]

# Used by type checkers for checks involving None (does not exist at runtime)
@final
class NoneType:
    def __bool__(self) -> Literal[False]: ...
