import _compression
import sys
from _compression import BaseStream
from _typeshed import ReadableBuffer, Self, StrOrBytesPath, WriteableBuffer
from typing import IO, Any, Iterable, Protocol, TextIO, TypeVar, overload
from typing_extensions import Literal, SupportsIndex, final

# The following attributes and methods are optional:
# def fileno(self) -> int: ...
# def close(self) -> object: ...
class _ReadableFileobj(_compression._Reader, Protocol): ...

class _WritableFileobj(Protocol):
    def write(self, __b: bytes) -> object: ...
    # The following attributes and methods are optional:
    # def fileno(self) -> int: ...
    # def close(self) -> object: ...

_T = TypeVar("_T")

def compress(data: bytes, compresslevel: int = ...) -> bytes: ...
def decompress(data: bytes) -> bytes: ...

_ReadBinaryMode = Literal["", "r", "rb"]
_WriteBinaryMode = Literal["w", "wb", "x", "xb", "a", "ab"]
_ReadTextMode = Literal["rt"]
_WriteTextMode = Literal["wt", "xt", "at"]

@overload
def open(
    filename: _ReadableFileobj,
    mode: _ReadBinaryMode = ...,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> BZ2File: ...
@overload
def open(
    filename: _ReadableFileobj,
    mode: _ReadTextMode,
    compresslevel: int = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...
@overload
def open(
    filename: _WritableFileobj,
    mode: _WriteBinaryMode,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> BZ2File: ...
@overload
def open(
    filename: _WritableFileobj,
    mode: _WriteTextMode,
    compresslevel: int = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...
@overload
def open(
    filename: StrOrBytesPath,
    mode: _ReadBinaryMode | _WriteBinaryMode = ...,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> BZ2File: ...
@overload
def open(
    filename: StrOrBytesPath,
    mode: _ReadTextMode | _WriteTextMode,
    compresslevel: int = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...

class BZ2File(BaseStream, IO[bytes]):
    def __enter__(self: Self) -> Self: ...
    if sys.version_info >= (3, 9):
        @overload
        def __init__(self, filename: _WritableFileobj, mode: _WriteBinaryMode, *, compresslevel: int = ...) -> None: ...
        @overload
        def __init__(self, filename: _ReadableFileobj, mode: _ReadBinaryMode = ..., *, compresslevel: int = ...) -> None: ...
        @overload
        def __init__(
            self, filename: StrOrBytesPath, mode: _ReadBinaryMode | _WriteBinaryMode = ..., *, compresslevel: int = ...
        ) -> None: ...
    else:
        @overload
        def __init__(
            self, filename: _WritableFileobj, mode: _WriteBinaryMode, buffering: Any | None = ..., compresslevel: int = ...
        ) -> None: ...
        @overload
        def __init__(
            self, filename: _ReadableFileobj, mode: _ReadBinaryMode = ..., buffering: Any | None = ..., compresslevel: int = ...
        ) -> None: ...
        @overload
        def __init__(
            self,
            filename: StrOrBytesPath,
            mode: _ReadBinaryMode | _WriteBinaryMode = ...,
            buffering: Any | None = ...,
            compresslevel: int = ...,
        ) -> None: ...
    def read(self, size: int | None = ...) -> bytes: ...
    def read1(self, size: int = ...) -> bytes: ...
    def readline(self, size: SupportsIndex = ...) -> bytes: ...  # type: ignore
    def readinto(self, b: WriteableBuffer) -> int: ...
    def readlines(self, size: SupportsIndex = ...) -> list[bytes]: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def write(self, data: ReadableBuffer) -> int: ...
    def writelines(self, seq: Iterable[ReadableBuffer]) -> None: ...

@final
class BZ2Compressor(object):
    def __init__(self, compresslevel: int = ...) -> None: ...
    def compress(self, __data: bytes) -> bytes: ...
    def flush(self) -> bytes: ...

@final
class BZ2Decompressor(object):
    def decompress(self, data: bytes, max_length: int = ...) -> bytes: ...
    @property
    def eof(self) -> bool: ...
    @property
    def needs_input(self) -> bool: ...
    @property
    def unused_data(self) -> bytes: ...
