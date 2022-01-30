import _compression
import sys
import zlib
from _typeshed import ReadableBuffer, StrOrBytesPath
from io import FileIO
from typing import Any, Protocol, TextIO, overload
from typing_extensions import Literal

_ReadBinaryMode = Literal["r", "rb"]
_WriteBinaryMode = Literal["a", "ab", "w", "wb", "x", "xb"]
_OpenTextMode = Literal["rt", "at", "wt", "xt"]

READ: Literal[1]
WRITE: Literal[2]

class _ReadableFileobj(Protocol):
    def read(self, __n: int) -> bytes: ...
    def seek(self, __n: int) -> Any: ...
    # The following attributes and methods are optional:
    # name: str
    # mode: str
    # def fileno() -> int: ...

class _WritableFileobj(Protocol):
    def write(self, __b: bytes) -> Any: ...
    def flush(self) -> Any: ...
    # The following attributes and methods are optional:
    # name: str
    # mode: str
    # def fileno() -> int: ...

@overload
def open(
    filename: StrOrBytesPath | _ReadableFileobj,
    mode: _ReadBinaryMode = ...,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> GzipFile: ...
@overload
def open(
    filename: StrOrBytesPath | _WritableFileobj,
    mode: _WriteBinaryMode,
    compresslevel: int = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> GzipFile: ...
@overload
def open(
    filename: StrOrBytesPath,
    mode: _OpenTextMode,
    compresslevel: int = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...
@overload
def open(
    filename: StrOrBytesPath | _ReadableFileobj | _WritableFileobj,
    mode: str,
    compresslevel: int = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> GzipFile | TextIO: ...

class _PaddedFile:
    file: _ReadableFileobj
    def __init__(self, f: _ReadableFileobj, prepend: bytes = ...) -> None: ...
    def read(self, size: int) -> bytes: ...
    def prepend(self, prepend: bytes = ...) -> None: ...
    def seek(self, off: int) -> int: ...
    def seekable(self) -> bool: ...

if sys.version_info >= (3, 8):
    class BadGzipFile(OSError): ...

class GzipFile(_compression.BaseStream):
    myfileobj: FileIO | None
    mode: Literal[1, 2]
    name: str
    compress: zlib._Compress
    fileobj: _ReadableFileobj | _WritableFileobj
    @overload
    def __init__(
        self,
        filename: StrOrBytesPath | None,
        mode: _ReadBinaryMode,
        compresslevel: int = ...,
        fileobj: _ReadableFileobj | None = ...,
        mtime: float | None = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        mode: _ReadBinaryMode,
        compresslevel: int = ...,
        fileobj: _ReadableFileobj | None = ...,
        mtime: float | None = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        filename: StrOrBytesPath | None,
        mode: _WriteBinaryMode,
        compresslevel: int = ...,
        fileobj: _WritableFileobj | None = ...,
        mtime: float | None = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        mode: _WriteBinaryMode,
        compresslevel: int = ...,
        fileobj: _WritableFileobj | None = ...,
        mtime: float | None = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        filename: StrOrBytesPath | None = ...,
        mode: str | None = ...,
        compresslevel: int = ...,
        fileobj: _ReadableFileobj | _WritableFileobj | None = ...,
        mtime: float | None = ...,
    ) -> None: ...
    @property
    def filename(self) -> str: ...
    @property
    def mtime(self) -> int | None: ...
    crc: int
    def write(self, data: ReadableBuffer) -> int: ...
    def read(self, size: int | None = ...) -> bytes: ...
    def read1(self, size: int = ...) -> bytes: ...
    def peek(self, n: int) -> bytes: ...
    @property
    def closed(self) -> bool: ...
    def close(self) -> None: ...
    def flush(self, zlib_mode: int = ...) -> None: ...
    def fileno(self) -> int: ...
    def rewind(self) -> None: ...
    def readable(self) -> bool: ...
    def writable(self) -> bool: ...
    def seekable(self) -> bool: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def readline(self, size: int | None = ...) -> bytes: ...

class _GzipReader(_compression.DecompressReader):
    def __init__(self, fp: _ReadableFileobj) -> None: ...
    def read(self, size: int = ...) -> bytes: ...

if sys.version_info >= (3, 8):
    def compress(data: bytes, compresslevel: int = ..., *, mtime: float | None = ...) -> bytes: ...

else:
    def compress(data: bytes, compresslevel: int = ...) -> bytes: ...

def decompress(data: bytes) -> bytes: ...
