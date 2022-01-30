import io
from types import TracebackType
from typing import IO, Callable, Dict, Iterable, Iterator, List, Mapping, Text, Tuple, Type

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
PAX_NUMBER_FIELDS: Dict[str, type]

ENCODING: str

TAR_PLAIN: int
TAR_GZIPPED: int

def open(
    name: Text | None = ...,
    mode: str = ...,
    fileobj: IO[bytes] | None = ...,
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

class TarFile(Iterable[TarInfo]):
    OPEN_METH: Mapping[str, str]
    name: Text | None
    mode: str
    fileobj: IO[bytes] | None
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
    posix: bool
    def __init__(
        self,
        name: Text | None = ...,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,
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
    def __enter__(self) -> TarFile: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    def __iter__(self) -> Iterator[TarInfo]: ...
    @classmethod
    def open(
        cls,
        name: Text | None = ...,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,
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
    ) -> TarFile: ...
    @classmethod
    def taropen(
        cls,
        name: Text | None,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,
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
    ) -> TarFile: ...
    @classmethod
    def gzopen(
        cls,
        name: Text | None,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,
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
    ) -> TarFile: ...
    @classmethod
    def bz2open(
        cls,
        name: Text | None,
        mode: str = ...,
        fileobj: IO[bytes] | None = ...,
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
    ) -> TarFile: ...
    @classmethod
    def xzopen(
        cls,
        name: Text | None,
        mode: str = ...,
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
    ) -> TarFile: ...
    def getmember(self, name: str) -> TarInfo: ...
    def getmembers(self) -> List[TarInfo]: ...
    def getnames(self) -> List[str]: ...
    def list(self, verbose: bool = ...) -> None: ...
    def next(self) -> TarInfo | None: ...
    def extractall(self, path: Text = ..., members: Iterable[TarInfo] | None = ...) -> None: ...
    def extract(self, member: str | TarInfo, path: Text = ...) -> None: ...
    def extractfile(self, member: str | TarInfo) -> IO[bytes] | None: ...
    def makedir(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def makefile(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def makeunknown(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def makefifo(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def makedev(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def makelink(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def chown(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def chmod(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def utime(self, tarinfo: TarInfo, targetpath: Text) -> None: ...  # undocumented
    def add(
        self,
        name: str,
        arcname: str | None = ...,
        recursive: bool = ...,
        exclude: Callable[[str], bool] | None = ...,
        filter: Callable[[TarInfo], TarInfo | None] | None = ...,
    ) -> None: ...
    def addfile(self, tarinfo: TarInfo, fileobj: IO[bytes] | None = ...) -> None: ...
    def gettarinfo(self, name: str | None = ..., arcname: str | None = ..., fileobj: IO[bytes] | None = ...) -> TarInfo: ...
    def close(self) -> None: ...

def is_tarfile(name: Text) -> bool: ...
def filemode(mode: int) -> str: ...  # undocumented

class TarFileCompat:
    def __init__(self, filename: str, mode: str = ..., compression: int = ...) -> None: ...

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
    def frombuf(cls, buf: bytes) -> TarInfo: ...
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
