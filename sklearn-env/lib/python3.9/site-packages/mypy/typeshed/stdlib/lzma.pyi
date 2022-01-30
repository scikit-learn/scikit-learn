import io
from _typeshed import ReadableBuffer, Self, StrOrBytesPath
from typing import IO, Any, Mapping, Sequence, TextIO, Union, overload
from typing_extensions import Literal, final

_OpenBinaryWritingMode = Literal["w", "wb", "x", "xb", "a", "ab"]
_OpenTextWritingMode = Literal["wt", "xt", "at"]

_PathOrFile = Union[StrOrBytesPath, IO[bytes]]

_FilterChain = Sequence[Mapping[str, Any]]

FORMAT_AUTO: int
FORMAT_XZ: int
FORMAT_ALONE: int
FORMAT_RAW: int
CHECK_NONE: int
CHECK_CRC32: int
CHECK_CRC64: int
CHECK_SHA256: int
CHECK_ID_MAX: int
CHECK_UNKNOWN: int
FILTER_LZMA1: int
FILTER_LZMA2: int
FILTER_DELTA: int
FILTER_X86: int
FILTER_IA64: int
FILTER_ARM: int
FILTER_ARMTHUMB: int
FILTER_SPARC: int
FILTER_POWERPC: int
MF_HC3: int
MF_HC4: int
MF_BT2: int
MF_BT3: int
MF_BT4: int
MODE_FAST: int
MODE_NORMAL: int
PRESET_DEFAULT: int
PRESET_EXTREME: int

# from _lzma.c
@final
class LZMADecompressor(object):
    def __init__(self, format: int | None = ..., memlimit: int | None = ..., filters: _FilterChain | None = ...) -> None: ...
    def decompress(self, data: bytes, max_length: int = ...) -> bytes: ...
    @property
    def check(self) -> int: ...
    @property
    def eof(self) -> bool: ...
    @property
    def unused_data(self) -> bytes: ...
    @property
    def needs_input(self) -> bool: ...

# from _lzma.c
@final
class LZMACompressor(object):
    def __init__(
        self, format: int | None = ..., check: int = ..., preset: int | None = ..., filters: _FilterChain | None = ...
    ) -> None: ...
    def compress(self, __data: bytes) -> bytes: ...
    def flush(self) -> bytes: ...

class LZMAError(Exception): ...

class LZMAFile(io.BufferedIOBase, IO[bytes]):
    def __init__(
        self,
        filename: _PathOrFile | None = ...,
        mode: str = ...,
        *,
        format: int | None = ...,
        check: int = ...,
        preset: int | None = ...,
        filters: _FilterChain | None = ...,
    ) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def close(self) -> None: ...
    @property
    def closed(self) -> bool: ...
    def fileno(self) -> int: ...
    def seekable(self) -> bool: ...
    def readable(self) -> bool: ...
    def writable(self) -> bool: ...
    def peek(self, size: int = ...) -> bytes: ...
    def read(self, size: int | None = ...) -> bytes: ...
    def read1(self, size: int = ...) -> bytes: ...
    def readline(self, size: int | None = ...) -> bytes: ...
    def write(self, data: ReadableBuffer) -> int: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def tell(self) -> int: ...

@overload
def open(
    filename: _PathOrFile,
    mode: Literal["r", "rb"] = ...,
    *,
    format: int | None = ...,
    check: Literal[-1] = ...,
    preset: None = ...,
    filters: _FilterChain | None = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> LZMAFile: ...
@overload
def open(
    filename: _PathOrFile,
    mode: _OpenBinaryWritingMode,
    *,
    format: int | None = ...,
    check: int = ...,
    preset: int | None = ...,
    filters: _FilterChain | None = ...,
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> LZMAFile: ...
@overload
def open(
    filename: StrOrBytesPath,
    mode: Literal["rt"],
    *,
    format: int | None = ...,
    check: Literal[-1] = ...,
    preset: None = ...,
    filters: _FilterChain | None = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...
@overload
def open(
    filename: StrOrBytesPath,
    mode: _OpenTextWritingMode,
    *,
    format: int | None = ...,
    check: int = ...,
    preset: int | None = ...,
    filters: _FilterChain | None = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> TextIO: ...
@overload
def open(
    filename: _PathOrFile,
    mode: str,
    *,
    format: int | None = ...,
    check: int = ...,
    preset: int | None = ...,
    filters: _FilterChain | None = ...,
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> LZMAFile | TextIO: ...
def compress(
    data: bytes, format: int = ..., check: int = ..., preset: int | None = ..., filters: _FilterChain | None = ...
) -> bytes: ...
def decompress(data: bytes, format: int = ..., memlimit: int | None = ..., filters: _FilterChain | None = ...) -> bytes: ...
def is_check_supported(__check_id: int) -> bool: ...
