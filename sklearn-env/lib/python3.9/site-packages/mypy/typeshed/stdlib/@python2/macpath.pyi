from genericpath import (
    commonprefix as commonprefix,
    exists as exists,
    getatime as getatime,
    getctime as getctime,
    getmtime as getmtime,
    getsize as getsize,
    isdir as isdir,
    isfile as isfile,
)

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
from typing import AnyStr, Text, Tuple, overload

altsep: str | None

def basename(s: AnyStr) -> AnyStr: ...
def dirname(s: AnyStr) -> AnyStr: ...
def normcase(path: AnyStr) -> AnyStr: ...
def normpath(s: AnyStr) -> AnyStr: ...
def realpath(path: AnyStr) -> AnyStr: ...
def islink(s: Text) -> bool: ...

# Make sure signatures are disjunct, and allow combinations of bytes and unicode.
# (Since Python 2 allows that, too)
# Note that e.g. os.path.join("a", "b", "c", "d", u"e") will still result in
# a type error.
@overload
def join(__p1: bytes, *p: bytes) -> bytes: ...
@overload
def join(__p1: bytes, __p2: bytes, __p3: bytes, __p4: Text, *p: Text) -> Text: ...
@overload
def join(__p1: bytes, __p2: bytes, __p3: Text, *p: Text) -> Text: ...
@overload
def join(__p1: bytes, __p2: Text, *p: Text) -> Text: ...
@overload
def join(__p1: Text, *p: Text) -> Text: ...
def split(s: AnyStr) -> Tuple[AnyStr, AnyStr]: ...
