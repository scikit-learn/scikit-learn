from _typeshed import BytesPath, StrOrBytesPath, StrPath
from genericpath import (
    commonprefix as commonprefix,
    exists as exists,
    getatime as getatime,
    getctime as getctime,
    getmtime as getmtime,
    getsize as getsize,
    isdir as isdir,
    isfile as isfile,
    samefile as samefile,
    sameopenfile as sameopenfile,
    samestat as samestat,
)
from os import PathLike

# Re-export common definitions from posixpath to reduce duplication
from posixpath import (
    abspath as abspath,
    curdir as curdir,
    defpath as defpath,
    devnull as devnull,
    expanduser as expanduser,
    expandvars as expandvars,
    extsep as extsep,
    isabs as isabs,
    lexists as lexists,
    pardir as pardir,
    pathsep as pathsep,
    sep as sep,
    splitdrive as splitdrive,
    splitext as splitext,
    supports_unicode_filenames as supports_unicode_filenames,
)
from typing import AnyStr, overload

altsep: str | None

@overload
def basename(s: PathLike[AnyStr]) -> AnyStr: ...
@overload
def basename(s: AnyStr) -> AnyStr: ...
@overload
def dirname(s: PathLike[AnyStr]) -> AnyStr: ...
@overload
def dirname(s: AnyStr) -> AnyStr: ...
@overload
def normcase(path: PathLike[AnyStr]) -> AnyStr: ...
@overload
def normcase(path: AnyStr) -> AnyStr: ...
@overload
def normpath(s: PathLike[AnyStr]) -> AnyStr: ...
@overload
def normpath(s: AnyStr) -> AnyStr: ...
@overload
def realpath(path: PathLike[AnyStr]) -> AnyStr: ...
@overload
def realpath(path: AnyStr) -> AnyStr: ...
def islink(s: StrOrBytesPath) -> bool: ...

# Mypy complains that the signatures overlap, but things seem to behave correctly anyway.
@overload
def join(s: StrPath, *paths: StrPath) -> str: ...
@overload
def join(s: BytesPath, *paths: BytesPath) -> bytes: ...
@overload
def split(s: PathLike[AnyStr]) -> tuple[AnyStr, AnyStr]: ...
@overload
def split(s: AnyStr) -> tuple[AnyStr, AnyStr]: ...
