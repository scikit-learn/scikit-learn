import io
import sys
from _typeshed import Self, StrPath
from os import PathLike
from types import TracebackType
from typing import IO, Any, Callable, Iterable, Iterator, Protocol, Sequence, Tuple, Type, overload
from typing_extensions import Literal

_DateTuple = Tuple[int, int, int, int, int, int]
_ReadWriteMode = Literal["r", "w"]
_ReadWriteBinaryMode = Literal["r", "w", "rb", "wb"]
_ZipFileMode = Literal["r", "w", "x", "a"]

class BadZipFile(Exception): ...

BadZipfile = BadZipFile
error = BadZipfile

class LargeZipFile(Exception): ...

class _ZipStream(Protocol):
    def read(self, __n: int) -> bytes: ...
    # The following methods are optional:
    # def seekable(self) -> bool: ...
    # def tell(self) -> int: ...
    # def seek(self, __n: int) -> object: ...

class _ClosableZipStream(_ZipStream, Protocol):
    def close(self) -> object: ...

class ZipExtFile(io.BufferedIOBase):
    MAX_N: int
    MIN_READ_SIZE: int

    if sys.version_info >= (3, 7):
        MAX_SEEK_READ: int

    newlines: list[bytes] | None
    mode: _ReadWriteMode
    name: str
    if sys.version_info >= (3, 7):
        @overload
        def __init__(
            self,
            fileobj: _ClosableZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            pwd: bytes | None,
            close_fileobj: Literal[True],
        ) -> None: ...
        @overload
        def __init__(
            self,
            fileobj: _ClosableZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            pwd: bytes | None = ...,
            *,
            close_fileobj: Literal[True],
        ) -> None: ...
        @overload
        def __init__(
            self,
            fileobj: _ZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            pwd: bytes | None = ...,
            close_fileobj: Literal[False] = ...,
        ) -> None: ...
    else:
        @overload
        def __init__(
            self,
            fileobj: _ClosableZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            decrypter: Callable[[Sequence[int]], bytes] | None,
            close_fileobj: Literal[True],
        ) -> None: ...
        @overload
        def __init__(
            self,
            fileobj: _ClosableZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            decrypter: Callable[[Sequence[int]], bytes] | None = ...,
            *,
            close_fileobj: Literal[True],
        ) -> None: ...
        @overload
        def __init__(
            self,
            fileobj: _ZipStream,
            mode: _ReadWriteMode,
            zipinfo: ZipInfo,
            decrypter: Callable[[Sequence[int]], bytes] | None = ...,
            close_fileobj: Literal[False] = ...,
        ) -> None: ...
    def read(self, n: int | None = ...) -> bytes: ...
    def readline(self, limit: int = ...) -> bytes: ...  # type: ignore
    def __repr__(self) -> str: ...
    def peek(self, n: int = ...) -> bytes: ...
    def read1(self, n: int | None) -> bytes: ...  # type: ignore
    if sys.version_info >= (3, 7):
        def seek(self, offset: int, whence: int = ...) -> int: ...

class _Writer(Protocol):
    def write(self, __s: str) -> object: ...

class ZipFile:
    filename: str | None
    debug: int
    comment: bytes
    filelist: list[ZipInfo]
    fp: IO[bytes] | None
    NameToInfo: dict[str, ZipInfo]
    start_dir: int  # undocumented
    compression: int  # undocumented
    compresslevel: int | None  # undocumented
    mode: _ZipFileMode  # undocumented
    pwd: str | None  # undocumented
    if sys.version_info >= (3, 8):
        def __init__(
            self,
            file: StrPath | IO[bytes],
            mode: _ZipFileMode = ...,
            compression: int = ...,
            allowZip64: bool = ...,
            compresslevel: int | None = ...,
            *,
            strict_timestamps: bool = ...,
        ) -> None: ...
    elif sys.version_info >= (3, 7):
        def __init__(
            self,
            file: StrPath | IO[bytes],
            mode: _ZipFileMode = ...,
            compression: int = ...,
            allowZip64: bool = ...,
            compresslevel: int | None = ...,
        ) -> None: ...
    else:
        def __init__(
            self, file: StrPath | IO[bytes], mode: _ZipFileMode = ..., compression: int = ..., allowZip64: bool = ...
        ) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    def close(self) -> None: ...
    def getinfo(self, name: str) -> ZipInfo: ...
    def infolist(self) -> list[ZipInfo]: ...
    def namelist(self) -> list[str]: ...
    def open(
        self, name: str | ZipInfo, mode: _ReadWriteMode = ..., pwd: bytes | None = ..., *, force_zip64: bool = ...
    ) -> IO[bytes]: ...
    def extract(self, member: str | ZipInfo, path: StrPath | None = ..., pwd: bytes | None = ...) -> str: ...
    def extractall(
        self, path: StrPath | None = ..., members: Iterable[str | ZipInfo] | None = ..., pwd: bytes | None = ...
    ) -> None: ...
    def printdir(self, file: _Writer | None = ...) -> None: ...
    def setpassword(self, pwd: bytes) -> None: ...
    def read(self, name: str | ZipInfo, pwd: bytes | None = ...) -> bytes: ...
    def testzip(self) -> str | None: ...
    if sys.version_info >= (3, 7):
        def write(
            self,
            filename: StrPath,
            arcname: StrPath | None = ...,
            compress_type: int | None = ...,
            compresslevel: int | None = ...,
        ) -> None: ...
    else:
        def write(self, filename: StrPath, arcname: StrPath | None = ..., compress_type: int | None = ...) -> None: ...
    if sys.version_info >= (3, 7):
        def writestr(
            self,
            zinfo_or_arcname: str | ZipInfo,
            data: bytes | str,
            compress_type: int | None = ...,
            compresslevel: int | None = ...,
        ) -> None: ...
    else:
        def writestr(self, zinfo_or_arcname: str | ZipInfo, data: bytes | str, compress_type: int | None = ...) -> None: ...

class PyZipFile(ZipFile):
    def __init__(
        self, file: str | IO[bytes], mode: _ZipFileMode = ..., compression: int = ..., allowZip64: bool = ..., optimize: int = ...
    ) -> None: ...
    def writepy(self, pathname: str, basename: str = ..., filterfunc: Callable[[str], bool] | None = ...) -> None: ...

class ZipInfo:
    filename: str
    date_time: _DateTuple
    compress_type: int
    comment: bytes
    extra: bytes
    create_system: int
    create_version: int
    extract_version: int
    reserved: int
    flag_bits: int
    volume: int
    internal_attr: int
    external_attr: int
    header_offset: int
    CRC: int
    compress_size: int
    file_size: int
    orig_filename: str  # undocumented
    def __init__(self, filename: str = ..., date_time: _DateTuple = ...) -> None: ...
    if sys.version_info >= (3, 8):
        @classmethod
        def from_file(cls, filename: StrPath, arcname: StrPath | None = ..., *, strict_timestamps: bool = ...) -> ZipInfo: ...
    else:
        @classmethod
        def from_file(cls, filename: StrPath, arcname: StrPath | None = ...) -> ZipInfo: ...
    def is_dir(self) -> bool: ...
    def FileHeader(self, zip64: bool | None = ...) -> bytes: ...

class _PathOpenProtocol(Protocol):
    def __call__(self, mode: _ReadWriteMode = ..., pwd: bytes | None = ..., *, force_zip64: bool = ...) -> IO[bytes]: ...

if sys.version_info >= (3, 8):
    class Path:
        @property
        def name(self) -> str: ...
        @property
        def parent(self) -> PathLike[str]: ...  # undocumented
        if sys.version_info >= (3, 10):
            @property
            def filename(self) -> PathLike[str]: ...  # undocumented
        def __init__(self, root: ZipFile | StrPath | IO[bytes], at: str = ...) -> None: ...
        if sys.version_info >= (3, 9):
            def open(self, mode: _ReadWriteBinaryMode = ..., *args: Any, pwd: bytes | None = ..., **kwargs: Any) -> IO[bytes]: ...
        else:
            @property
            def open(self) -> _PathOpenProtocol: ...
        def iterdir(self) -> Iterator[Path]: ...
        def is_dir(self) -> bool: ...
        def is_file(self) -> bool: ...
        def exists(self) -> bool: ...
        def read_text(
            self,
            encoding: str | None = ...,
            errors: str | None = ...,
            newline: str | None = ...,
            line_buffering: bool = ...,
            write_through: bool = ...,
        ) -> str: ...
        def read_bytes(self) -> bytes: ...
        if sys.version_info >= (3, 10):
            def joinpath(self, *other: StrPath) -> Path: ...
        else:
            def joinpath(self, add: StrPath) -> Path: ...  # undocumented
        def __truediv__(self, add: StrPath) -> Path: ...

def is_zipfile(filename: StrPath | IO[bytes]) -> bool: ...

ZIP_STORED: int
ZIP_DEFLATED: int
ZIP64_LIMIT: int
ZIP_FILECOUNT_LIMIT: int
ZIP_MAX_COMMENT: int
ZIP_BZIP2: int
ZIP_LZMA: int
