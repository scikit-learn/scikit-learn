import zipfile
import types
from _typeshed import StrOrBytesPath, StrPath, SupportsRead, SupportsWrite, SupportsKeysAndGetItem
from re import Pattern
from collections.abc import Collection, Mapping, Iterator, Sequence, Callable, Iterable
from typing import (
    Literal as L,
    Any,
    TypeVar,
    Generic,
    IO,
    overload,
    Protocol,
    type_check_only,
)
from typing_extensions import deprecated

from numpy import (
    recarray,
    dtype,
    generic,
    float64,
    void,
    record,
)
from numpy.ma.mrecords import MaskedRecords
from numpy._core.multiarray import packbits, unpackbits
from numpy._typing import (
    ArrayLike,
    DTypeLike,
    NDArray,
    _DTypeLike,
    _SupportsArrayFunc,
)

__all__ = [
    "savetxt",
    "loadtxt",
    "genfromtxt",
    "load",
    "save",
    "savez",
    "savez_compressed",
    "packbits",
    "unpackbits",
    "fromregex",
]

_T = TypeVar("_T")
_T_contra = TypeVar("_T_contra", contravariant=True)
_T_co = TypeVar("_T_co", covariant=True)
_SCT = TypeVar("_SCT", bound=generic)

@type_check_only
class _SupportsReadSeek(SupportsRead[_T_co], Protocol[_T_co]):
    def seek(self, offset: int, whence: int, /) -> object: ...

class BagObj(Generic[_T_co]):
    def __init__(self, obj: SupportsKeysAndGetItem[str, _T_co]) -> None: ...
    def __getattribute__(self, key: str) -> _T_co: ...
    def __dir__(self) -> list[str]: ...

class NpzFile(Mapping[str, NDArray[Any]]):
    zip: zipfile.ZipFile
    fid: None | IO[str]
    files: list[str]
    allow_pickle: bool
    pickle_kwargs: None | Mapping[str, Any]
    _MAX_REPR_ARRAY_COUNT: int
    # Represent `f` as a mutable property so we can access the type of `self`
    @property
    def f(self: _T) -> BagObj[_T]: ...
    @f.setter
    def f(self: _T, value: BagObj[_T]) -> None: ...
    def __init__(
        self,
        fid: IO[str],
        own_fid: bool = ...,
        allow_pickle: bool = ...,
        pickle_kwargs: None | Mapping[str, Any] = ...,
    ) -> None: ...
    def __enter__(self: _T) -> _T: ...
    def __exit__(
        self,
        exc_type: None | type[BaseException],
        exc_value: None | BaseException,
        traceback: None | types.TracebackType,
        /,
    ) -> None: ...
    def close(self) -> None: ...
    def __del__(self) -> None: ...
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...
    def __getitem__(self, key: str) -> NDArray[Any]: ...
    def __contains__(self, key: str) -> bool: ...
    def __repr__(self) -> str: ...

class DataSource:
    def __init__(self, destpath: StrPath | None = ...) -> None: ...
    def __del__(self) -> None: ...
    def abspath(self, path: str) -> str: ...
    def exists(self, path: str) -> bool: ...

    # Whether the file-object is opened in string or bytes mode (by default)
    # depends on the file-extension of `path`
    def open(
        self,
        path: str,
        mode: str = ...,
        encoding: None | str = ...,
        newline: None | str = ...,
    ) -> IO[Any]: ...

# NOTE: Returns a `NpzFile` if file is a zip file;
# returns an `ndarray`/`memmap` otherwise
def load(
    file: StrOrBytesPath | _SupportsReadSeek[bytes],
    mmap_mode: L[None, "r+", "r", "w+", "c"] = ...,
    allow_pickle: bool = ...,
    fix_imports: bool = ...,
    encoding: L["ASCII", "latin1", "bytes"] = ...,
) -> Any: ...

@overload
def save(
    file: StrPath | SupportsWrite[bytes],
    arr: ArrayLike,
    allow_pickle: bool = ...,
) -> None: ...
@overload
@deprecated("The 'fix_imports' flag is deprecated in NumPy 2.1.")
def save(
    file: StrPath | SupportsWrite[bytes],
    arr: ArrayLike,
    allow_pickle: bool = ...,
    *,
    fix_imports: bool,
) -> None: ...
@overload
@deprecated("The 'fix_imports' flag is deprecated in NumPy 2.1.")
def save(
    file: StrPath | SupportsWrite[bytes],
    arr: ArrayLike,
    allow_pickle: bool,
    fix_imports: bool,
) -> None: ...

def savez(
    file: StrPath | SupportsWrite[bytes],
    *args: ArrayLike,
    allow_pickle: bool = ...,
    **kwds: ArrayLike,
) -> None: ...

def savez_compressed(
    file: StrPath | SupportsWrite[bytes],
    *args: ArrayLike,
    allow_pickle: bool = ...,
    **kwds: ArrayLike,
) -> None: ...

# File-like objects only have to implement `__iter__` and,
# optionally, `encoding`
@overload
def loadtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: None = ...,
    comments: None | str | Sequence[str] = ...,
    delimiter: None | str = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] | Callable[[str], Any] = ...,
    skiprows: int = ...,
    usecols: int | Sequence[int] | None = ...,
    unpack: bool = ...,
    ndmin: L[0, 1, 2] = ...,
    encoding: None | str = ...,
    max_rows: None | int = ...,
    *,
    quotechar: None | str = ...,
    like: None | _SupportsArrayFunc = ...
) -> NDArray[float64]: ...
@overload
def loadtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: _DTypeLike[_SCT],
    comments: None | str | Sequence[str] = ...,
    delimiter: None | str = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] | Callable[[str], Any]  = ...,
    skiprows: int = ...,
    usecols: int | Sequence[int] | None = ...,
    unpack: bool = ...,
    ndmin: L[0, 1, 2] = ...,
    encoding: None | str = ...,
    max_rows: None | int = ...,
    *,
    quotechar: None | str = ...,
    like: None | _SupportsArrayFunc = ...
) -> NDArray[_SCT]: ...
@overload
def loadtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: DTypeLike,
    comments: None | str | Sequence[str] = ...,
    delimiter: None | str = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] | Callable[[str], Any]  = ...,
    skiprows: int = ...,
    usecols: int | Sequence[int] | None = ...,
    unpack: bool = ...,
    ndmin: L[0, 1, 2] = ...,
    encoding: None | str = ...,
    max_rows: None | int = ...,
    *,
    quotechar: None | str = ...,
    like: None | _SupportsArrayFunc = ...
) -> NDArray[Any]: ...

def savetxt(
    fname: StrPath | SupportsWrite[str] | SupportsWrite[bytes],
    X: ArrayLike,
    fmt: str | Sequence[str] = ...,
    delimiter: str = ...,
    newline: str = ...,
    header: str = ...,
    footer: str = ...,
    comments: str = ...,
    encoding: None | str = ...,
) -> None: ...

@overload
def fromregex(
    file: StrPath | SupportsRead[str] | SupportsRead[bytes],
    regexp: str | bytes | Pattern[Any],
    dtype: _DTypeLike[_SCT],
    encoding: None | str = ...
) -> NDArray[_SCT]: ...
@overload
def fromregex(
    file: StrPath | SupportsRead[str] | SupportsRead[bytes],
    regexp: str | bytes | Pattern[Any],
    dtype: DTypeLike,
    encoding: None | str = ...
) -> NDArray[Any]: ...

@overload
def genfromtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: None = ...,
    comments: str = ...,
    delimiter: None | str | int | Iterable[int] = ...,
    skip_header: int = ...,
    skip_footer: int = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] = ...,
    missing_values: Any = ...,
    filling_values: Any = ...,
    usecols: None | Sequence[int] = ...,
    names: L[None, True] | str | Collection[str] = ...,
    excludelist: None | Sequence[str] = ...,
    deletechars: str = ...,
    replace_space: str = ...,
    autostrip: bool = ...,
    case_sensitive: bool | L['upper', 'lower'] = ...,
    defaultfmt: str = ...,
    unpack: None | bool = ...,
    usemask: bool = ...,
    loose: bool = ...,
    invalid_raise: bool = ...,
    max_rows: None | int = ...,
    encoding: str = ...,
    *,
    ndmin: L[0, 1, 2] = ...,
    like: None | _SupportsArrayFunc = ...,
) -> NDArray[Any]: ...
@overload
def genfromtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: _DTypeLike[_SCT],
    comments: str = ...,
    delimiter: None | str | int | Iterable[int] = ...,
    skip_header: int = ...,
    skip_footer: int = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] = ...,
    missing_values: Any = ...,
    filling_values: Any = ...,
    usecols: None | Sequence[int] = ...,
    names: L[None, True] | str | Collection[str] = ...,
    excludelist: None | Sequence[str] = ...,
    deletechars: str = ...,
    replace_space: str = ...,
    autostrip: bool = ...,
    case_sensitive: bool | L['upper', 'lower'] = ...,
    defaultfmt: str = ...,
    unpack: None | bool = ...,
    usemask: bool = ...,
    loose: bool = ...,
    invalid_raise: bool = ...,
    max_rows: None | int = ...,
    encoding: str = ...,
    *,
    ndmin: L[0, 1, 2] = ...,
    like: None | _SupportsArrayFunc = ...,
) -> NDArray[_SCT]: ...
@overload
def genfromtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    dtype: DTypeLike,
    comments: str = ...,
    delimiter: None | str | int | Iterable[int] = ...,
    skip_header: int = ...,
    skip_footer: int = ...,
    converters: None | Mapping[int | str, Callable[[str], Any]] = ...,
    missing_values: Any = ...,
    filling_values: Any = ...,
    usecols: None | Sequence[int] = ...,
    names: L[None, True] | str | Collection[str] = ...,
    excludelist: None | Sequence[str] = ...,
    deletechars: str = ...,
    replace_space: str = ...,
    autostrip: bool = ...,
    case_sensitive: bool | L['upper', 'lower'] = ...,
    defaultfmt: str = ...,
    unpack: None | bool = ...,
    usemask: bool = ...,
    loose: bool = ...,
    invalid_raise: bool = ...,
    max_rows: None | int = ...,
    encoding: str = ...,
    *,
    ndmin: L[0, 1, 2] = ...,
    like: None | _SupportsArrayFunc = ...,
) -> NDArray[Any]: ...

@overload
def recfromtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    *,
    usemask: L[False] = ...,
    **kwargs: Any,
) -> recarray[Any, dtype[record]]: ...
@overload
def recfromtxt(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    *,
    usemask: L[True],
    **kwargs: Any,
) -> MaskedRecords[Any, dtype[void]]: ...

@overload
def recfromcsv(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    *,
    usemask: L[False] = ...,
    **kwargs: Any,
) -> recarray[Any, dtype[record]]: ...
@overload
def recfromcsv(
    fname: StrPath | Iterable[str] | Iterable[bytes],
    *,
    usemask: L[True],
    **kwargs: Any,
) -> MaskedRecords[Any, dtype[void]]: ...
