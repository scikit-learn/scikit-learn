import bz2
import io
import sys
from _typeshed import Self, StrOrBytesPath, StrPath
from collections.abc import Callable, Iterable, Iterator, Mapping
from gzip import _ReadableFileobj as _GzipReadableFileobj, _WritableFileobj as _GzipWritableFileobj
from types import TracebackType
from typing import IO, Protocol, Tuple, Type, TypeVar, overload
from typing_extensions import Literal

_TF = TypeVar("_TF", bound=TarFile)

class _Fileobj(Protocol):
    def read(self, __size: int) -> bytes: ...
    def write(self, __b: bytes) -> object: ...
    def tell(self) -> int: ...
    def seek(self, __pos: int) -> object: ...
    def close(self) -> object: ...
    # Optional fields:
    # name: str | bytes
    # mode: Literal["rb", "r+b", "wb", "xb"]

class _Bz2ReadableFileobj(bz2._ReadableFileobj):
    def close(self) -> object: ...

class _Bz2WritableFileobj(bz2._WritableFileobj):
    def close(self) -> object: ...

# tar constants
NUL: bytes
BLOCKSIZE: int
RECORDSIZE: int
GNU_MAGIC: bytes
POSIX_MAGIC: bytes

LENGTH_NAME: int
LENGTH_LINK: int
LENGTH_PREFIX: int

REGTYPE: bytes
AREGTYPE: bytes
LNKTYPE: bytes
SYMTYPE: bytes
CONTTYPE: bytes
BLKTYPE: bytes
DIRTYPE: bytes
FIFOTYPE: bytes
CHRTYPE: bytes

GNUTYPE_LONGNAME: bytes
GNUTYPE_LONGLINK: bytes
GNUTYPE_SPARSE: bytes

XHDTYPE: bytes
XGLTYPE: bytes
SOLARIS_XHDTYPE: bytes

USTAR_FORMAT: int
GNU_FORMAT: int
PAX_FORMAT: int
DEFAULT_FORMAT: int

# tarfile constants

SUPPORTED_TYPES: Tuple[bytes, ...]
REGULAR_TYPES: Tuple[bytes, ...]
GNU_TYPES: Tuple[bytes, ...]
PAX_FIELDS: Tuple[str, ...]
PAX_NUMBER_FIELDS: dict[str, type]
PAX_NAME_FIELDS: set[str]

ENCODING: str

def open(
    name: StrOrBytesPath | None = ...,
    mode: str = ...,
    fileobj: IO[bytes] | None = ...,  # depends on mode
    bufsize: int = ...,
    *,
    format: int | None = ...,
    tarinfo: Type[TarInfo] | None = ...,
    dereference: bool | None = ...,
    ignore_zeros: bool | None = ...,
    encoding: str | None = ...,
    errors: str = ...,
    pax_headers: Mapping[str, str] | None = ...,
    debug: int | None = ...,
    errorlevel: int | None = ...,
    compresslevel: int | None = ...,
) -> TarFile: ...

class ExFileObject(io.BufferedReader):
    def __init__(self, tarfile: TarFile, tarinfo: TarInfo) -> None: ...

_list = list  # conflicts with method name

class TarFile:
    OPEN_METH: Mapping[str, str]
    name: StrOrBytesPath | None
    mode: Literal["r", "a", "w", "x"]
    fileobj: _Fileobj | None
    format: int | None
    tarinfo: Type[TarInfo]
    dereference: bool | None
    ignore_zeros: bool | None
    encoding: str | None
    errors: str
    fileobject: Type[ExFileObject]
    pax_headers: Mapping[str, str] | None
    debug: int | None
    errorlevel: int | None
    offset: int  # undocumented
    def __init__(
        self,
        name: StrOrBytesPath | None = ...,
        mode: Literal["r", "a", "w", "x"] = ...,
        fileobj: _Fileobj | None = ...,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        errors: str = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
        copybufsize: int | None = ...,  # undocumented
    ) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    def __iter__(self) -> Iterator[TarInfo]: ...
    @classmethod
    def open(
        cls: Type[_TF],
        name: StrOrBytesPath | None = ...,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,  # depends on mode
        bufsize: int = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        errors: str = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @classmethod
    def taropen(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["r", "a", "w", "x"] = ...,
        fileobj: _Fileobj | None = ...,
        *,
        compresslevel: int = ...,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @overload
    @classmethod
    def gzopen(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["r"] = ...,
        fileobj: _GzipReadableFileobj | None = ...,
        compresslevel: int = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @overload
    @classmethod
    def gzopen(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["w", "x"],
        fileobj: _GzipWritableFileobj | None = ...,
        compresslevel: int = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @overload
    @classmethod
    def bz2open(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["w", "x"],
        fileobj: _Bz2WritableFileobj | None = ...,
        compresslevel: int = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @overload
    @classmethod
    def bz2open(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["r"] = ...,
        fileobj: _Bz2ReadableFileobj | None = ...,
        compresslevel: int = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    @classmethod
    def xzopen(
        cls: Type[_TF],
        name: StrOrBytesPath | None,
        mode: Literal["r", "w", "x"] = ...,
        fileobj: IO[bytes] | None = ...,
        preset: int | None = ...,
        *,
        format: int | None = ...,
        tarinfo: Type[TarInfo] | None = ...,
        dereference: bool | None = ...,
        ignore_zeros: bool | None = ...,
        encoding: str | None = ...,
        pax_headers: Mapping[str, str] | None = ...,
        debug: int | None = ...,
        errorlevel: int | None = ...,
    ) -> _TF: ...
    def getmember(self, name: str) -> TarInfo: ...
    def getmembers(self) -> _list[TarInfo]: ...
    def getnames(self) -> _list[str]: ...
    def list(self, verbose: bool = ..., *, members: _list[TarInfo] | None = ...) -> None: ...
    def next(self) -> TarInfo | None: ...
    def extractall(
        self, path: StrOrBytesPath = ..., members: Iterable[TarInfo] | None = ..., *, numeric_owner: bool = ...
    ) -> None: ...
    def extract(
        self, member: str | TarInfo, path: StrOrBytesPath = ..., set_attrs: bool = ..., *, numeric_owner: bool = ...
    ) -> None: ...
    def _extract_member(
        self, tarinfo: TarInfo, targetpath: str, set_attrs: bool = ..., numeric_owner: bool = ...
    ) -> None: ...  # undocumented
    def extractfile(self, member: str | TarInfo) -> IO[bytes] | None: ...
    def makedir(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def makefile(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def makeunknown(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def makefifo(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def makedev(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def makelink(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def chown(self, tarinfo: TarInfo, targetpath: StrOrBytesPath, numeric_owner: bool) -> None: ...  # undocumented
    def chmod(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    def utime(self, tarinfo: TarInfo, targetpath: StrOrBytesPath) -> None: ...  # undocumented
    if sys.version_info >= (3, 7):
        def add(
            self,
            name: StrPath,
            arcname: StrPath | None = ...,
            recursive: bool = ...,
            *,
            filter: Callable[[TarInfo], TarInfo | None] | None = ...,
        ) -> None: ...
    else:
        def add(
            self,
            name: StrPath,
            arcname: StrPath | None = ...,
            recursive: bool = ...,
            exclude: Callable[[str], bool] | None = ...,
            *,
            filter: Callable[[TarInfo], TarInfo | None] | None = ...,
        ) -> None: ...
    def addfile(self, tarinfo: TarInfo, fileobj: IO[bytes] | None = ...) -> None: ...
    def gettarinfo(self, name: str | None = ..., arcname: str | None = ..., fileobj: IO[bytes] | None = ...) -> TarInfo: ...
    def close(self) -> None: ...

if sys.version_info >= (3, 9):
    def is_tarfile(name: StrOrBytesPath | IO[bytes]) -> bool: ...

else:
    def is_tarfile(name: StrOrBytesPath) -> bool: ...

if sys.version_info < (3, 8):
    def filemode(mode: int) -> str: ...  # undocumented

class TarError(Exception): ...
class ReadError(TarError): ...
class CompressionError(TarError): ...
class StreamError(TarError): ...
class ExtractError(TarError): ...
class HeaderError(TarError): ...

class TarInfo:
    name: str
    path: str
    size: int
    mtime: int
    chksum: int
    devmajor: int
    devminor: int
    offset: int
    offset_data: int
    sparse: bytes | None
    tarfile: TarFile | None
    mode: int
    type: bytes
    linkname: str
    uid: int
    gid: int
    uname: str
    gname: str
    pax_headers: Mapping[str, str]
    def __init__(self, name: str = ...) -> None: ...
    @classmethod
    def frombuf(cls, buf: bytes, encoding: str, errors: str) -> TarInfo: ...
    @classmethod
    def fromtarfile(cls, tarfile: TarFile) -> TarInfo: ...
    @property
    def linkpath(self) -> str: ...
    @linkpath.setter
    def linkpath(self, linkname: str) -> None: ...
    def get_info(self) -> Mapping[str, str | int | bytes | Mapping[str, str]]: ...
    def tobuf(self, format: int | None = ..., encoding: str | None = ..., errors: str = ...) -> bytes: ...
    def create_ustar_header(
        self, info: Mapping[str, str | int | bytes | Mapping[str, str]], encoding: str, errors: str
    ) -> bytes: ...
    def create_gnu_header(
        self, info: Mapping[str, str | int | bytes | Mapping[str, str]], encoding: str, errors: str
    ) -> bytes: ...
    def create_pax_header(self, info: Mapping[str, str | int | bytes | Mapping[str, str]], encoding: str) -> bytes: ...
    @classmethod
    def create_pax_global_header(cls, pax_headers: Mapping[str, str]) -> bytes: ...
    def isfile(self) -> bool: ...
    def isreg(self) -> bool: ...
    def issparse(self) -> bool: ...
    def isdir(self) -> bool: ...
    def issym(self) -> bool: ...
    def islnk(self) -> bool: ...
    def ischr(self) -> bool: ...
    def isblk(self) -> bool: ...
    def isfifo(self) -> bool: ...
    def isdev(self) -> bool: ...
