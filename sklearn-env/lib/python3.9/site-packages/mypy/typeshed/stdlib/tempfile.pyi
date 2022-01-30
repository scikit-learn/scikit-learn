import os
import sys
from _typeshed import Self
from types import TracebackType
from typing import IO, Any, AnyStr, Generic, Iterable, Iterator, Tuple, Type, Union, overload
from typing_extensions import Literal

if sys.version_info >= (3, 9):
    from types import GenericAlias

# global variables
TMP_MAX: int
tempdir: str | None
template: str

_DirT = Union[AnyStr, os.PathLike[AnyStr]]

if sys.version_info >= (3, 8):
    @overload
    def NamedTemporaryFile(
        mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"],
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
        *,
        errors: str | None = ...,
    ) -> _TemporaryFileWrapper[str]: ...
    @overload
    def NamedTemporaryFile(
        mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
        *,
        errors: str | None = ...,
    ) -> _TemporaryFileWrapper[bytes]: ...
    @overload
    def NamedTemporaryFile(
        mode: str = ...,
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
        *,
        errors: str | None = ...,
    ) -> _TemporaryFileWrapper[Any]: ...

else:
    @overload
    def NamedTemporaryFile(
        mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"],
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
    ) -> _TemporaryFileWrapper[str]: ...
    @overload
    def NamedTemporaryFile(
        mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
    ) -> _TemporaryFileWrapper[bytes]: ...
    @overload
    def NamedTemporaryFile(
        mode: str = ...,
        buffering: int = ...,
        encoding: str | None = ...,
        newline: str | None = ...,
        suffix: AnyStr | None = ...,
        prefix: AnyStr | None = ...,
        dir: _DirT[AnyStr] | None = ...,
        delete: bool = ...,
    ) -> _TemporaryFileWrapper[Any]: ...

if sys.platform == "win32":
    TemporaryFile = NamedTemporaryFile
else:
    if sys.version_info >= (3, 8):
        @overload
        def TemporaryFile(
            mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"],
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
            *,
            errors: str | None = ...,
        ) -> IO[str]: ...
        @overload
        def TemporaryFile(
            mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
            *,
            errors: str | None = ...,
        ) -> IO[bytes]: ...
        @overload
        def TemporaryFile(
            mode: str = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
            *,
            errors: str | None = ...,
        ) -> IO[Any]: ...
    else:
        @overload
        def TemporaryFile(
            mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"],
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
        ) -> IO[str]: ...
        @overload
        def TemporaryFile(
            mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
        ) -> IO[bytes]: ...
        @overload
        def TemporaryFile(
            mode: str = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
        ) -> IO[Any]: ...

class _TemporaryFileWrapper(Generic[AnyStr], IO[AnyStr]):
    file: IO[AnyStr]  # io.TextIOWrapper, io.BufferedReader or io.BufferedWriter
    name: str
    delete: bool
    def __init__(self, file: IO[AnyStr], name: str, delete: bool = ...) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(self, exc: Type[BaseException] | None, value: BaseException | None, tb: TracebackType | None) -> bool | None: ...
    def __getattr__(self, name: str) -> Any: ...
    def close(self) -> None: ...
    # These methods don't exist directly on this object, but
    # are delegated to the underlying IO object through __getattr__.
    # We need to add them here so that this class is concrete.
    def __iter__(self) -> Iterator[AnyStr]: ...
    # FIXME: __next__ doesn't actually exist on this class and should be removed:
    #        see also https://github.com/python/typeshed/pull/5456#discussion_r633068648
    # >>> import tempfile
    # >>> ntf=tempfile.NamedTemporaryFile()
    # >>> next(ntf)
    # Traceback (most recent call last):
    #   File "<stdin>", line 1, in <module>
    # TypeError: '_TemporaryFileWrapper' object is not an iterator
    def __next__(self) -> AnyStr: ...
    def fileno(self) -> int: ...
    def flush(self) -> None: ...
    def isatty(self) -> bool: ...
    def read(self, n: int = ...) -> AnyStr: ...
    def readable(self) -> bool: ...
    def readline(self, limit: int = ...) -> AnyStr: ...
    def readlines(self, hint: int = ...) -> list[AnyStr]: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def seekable(self) -> bool: ...
    def tell(self) -> int: ...
    def truncate(self, size: int | None = ...) -> int: ...
    def writable(self) -> bool: ...
    def write(self, s: AnyStr) -> int: ...
    def writelines(self, lines: Iterable[AnyStr]) -> None: ...

# It does not actually derive from IO[AnyStr], but it does implement the
# protocol.
class SpooledTemporaryFile(IO[AnyStr]):
    @property
    def encoding(self) -> str: ...  # undocumented
    @property
    def newlines(self) -> str | Tuple[str, ...] | None: ...  # undocumented
    # bytes needs to go first, as default mode is to open as bytes
    if sys.version_info >= (3, 8):
        @overload
        def __init__(
            self: SpooledTemporaryFile[bytes],
            max_size: int = ...,
            mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
            *,
            errors: str | None = ...,
        ) -> None: ...
        @overload
        def __init__(
            self: SpooledTemporaryFile[str],
            max_size: int = ...,
            mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
            *,
            errors: str | None = ...,
        ) -> None: ...
        @overload
        def __init__(
            self,
            max_size: int = ...,
            mode: str = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
            *,
            errors: str | None = ...,
        ) -> None: ...
        @property
        def errors(self) -> str | None: ...
    else:
        @overload
        def __init__(
            self: SpooledTemporaryFile[bytes],
            max_size: int = ...,
            mode: Literal["rb", "wb", "ab", "xb", "r+b", "w+b", "a+b", "x+b"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
        ) -> None: ...
        @overload
        def __init__(
            self: SpooledTemporaryFile[str],
            max_size: int = ...,
            mode: Literal["r", "w", "a", "x", "r+", "w+", "a+", "x+", "rt", "wt", "at", "xt", "r+t", "w+t", "a+t", "x+t"] = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
        ) -> None: ...
        @overload
        def __init__(
            self,
            max_size: int = ...,
            mode: str = ...,
            buffering: int = ...,
            encoding: str | None = ...,
            newline: str | None = ...,
            suffix: str | None = ...,
            prefix: str | None = ...,
            dir: str | None = ...,
        ) -> None: ...
    def rollover(self) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> bool | None: ...
    # These methods are copied from the abstract methods of IO, because
    # SpooledTemporaryFile implements IO.
    # See also https://github.com/python/typeshed/pull/2452#issuecomment-420657918.
    def close(self) -> None: ...
    def fileno(self) -> int: ...
    def flush(self) -> None: ...
    def isatty(self) -> bool: ...
    def read(self, n: int = ...) -> AnyStr: ...
    def readline(self, limit: int = ...) -> AnyStr: ...
    def readlines(self, hint: int = ...) -> list[AnyStr]: ...
    def seek(self, offset: int, whence: int = ...) -> int: ...
    def tell(self) -> int: ...
    def truncate(self, size: int | None = ...) -> int: ...
    def write(self, s: AnyStr) -> int: ...
    def writelines(self, iterable: Iterable[AnyStr]) -> None: ...
    def __iter__(self) -> Iterator[AnyStr]: ...
    # Other than the following methods, which do not exist on SpooledTemporaryFile
    def readable(self) -> bool: ...
    def seekable(self) -> bool: ...
    def writable(self) -> bool: ...
    def __next__(self) -> AnyStr: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

class TemporaryDirectory(Generic[AnyStr]):
    name: AnyStr
    if sys.version_info >= (3, 10):
        @overload
        def __init__(
            self: TemporaryDirectory[str],
            suffix: None = ...,
            prefix: None = ...,
            dir: None = ...,
            ignore_cleanup_errors: bool = ...,
        ) -> None: ...
        @overload
        def __init__(
            self,
            suffix: AnyStr | None = ...,
            prefix: AnyStr | None = ...,
            dir: _DirT[AnyStr] | None = ...,
            ignore_cleanup_errors: bool = ...,
        ) -> None: ...
    else:
        @overload
        def __init__(self: TemporaryDirectory[str], suffix: None = ..., prefix: None = ..., dir: None = ...) -> None: ...
        @overload
        def __init__(self, suffix: AnyStr | None = ..., prefix: AnyStr | None = ..., dir: _DirT[AnyStr] | None = ...) -> None: ...
    def cleanup(self) -> None: ...
    def __enter__(self) -> AnyStr: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

@overload
def mkstemp() -> tuple[int, str]: ...
@overload
def mkstemp(
    suffix: AnyStr | None = ..., prefix: AnyStr | None = ..., dir: _DirT[AnyStr] | None = ..., text: bool = ...
) -> tuple[int, AnyStr]: ...
@overload
def mkdtemp() -> str: ...
@overload
def mkdtemp(suffix: AnyStr | None = ..., prefix: AnyStr | None = ..., dir: _DirT[AnyStr] | None = ...) -> AnyStr: ...
@overload
def mktemp() -> str: ...
@overload
def mktemp(suffix: AnyStr | None = ..., prefix: AnyStr | None = ..., dir: _DirT[AnyStr] | None = ...) -> AnyStr: ...
def gettempdirb() -> bytes: ...
def gettempprefixb() -> bytes: ...
def gettempdir() -> str: ...
def gettempprefix() -> str: ...
