from typing import Optional, Union, Callable, Any, Literal, TypedDict

from numpy import _SupportsWrite

_ErrKind = Literal["ignore", "warn", "raise", "call", "print", "log"]
_ErrFunc = Callable[[str, int], Any]

class _ErrDict(TypedDict):
    divide: _ErrKind
    over: _ErrKind
    under: _ErrKind
    invalid: _ErrKind

class _ErrDictOptional(TypedDict, total=False):
    all: Optional[_ErrKind]
    divide: Optional[_ErrKind]
    over: Optional[_ErrKind]
    under: Optional[_ErrKind]
    invalid: Optional[_ErrKind]

def seterr(
    all: Optional[_ErrKind] = ...,
    divide: Optional[_ErrKind] = ...,
    over: Optional[_ErrKind] = ...,
    under: Optional[_ErrKind] = ...,
    invalid: Optional[_ErrKind] = ...,
) -> _ErrDict: ...
def geterr() -> _ErrDict: ...
def setbufsize(size: int) -> int: ...
def getbufsize() -> int: ...
def seterrcall(
    func: Union[None, _ErrFunc, _SupportsWrite[str]]
) -> Union[None, _ErrFunc, _SupportsWrite[str]]: ...
def geterrcall() -> Union[None, _ErrFunc, _SupportsWrite[str]]: ...

# See `numpy/__init__.pyi` for the `errstate` class
