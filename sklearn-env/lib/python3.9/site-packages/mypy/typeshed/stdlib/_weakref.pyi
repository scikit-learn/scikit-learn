import sys
from typing import Any, Callable, Generic, TypeVar, overload
from typing_extensions import final

if sys.version_info >= (3, 9):
    from types import GenericAlias

_C = TypeVar("_C", bound=Callable[..., Any])
_T = TypeVar("_T")

@final
class CallableProxyType(Generic[_C]):  # "weakcallableproxy"
    def __getattr__(self, attr: str) -> Any: ...

@final
class ProxyType(Generic[_T]):  # "weakproxy"
    def __getattr__(self, attr: str) -> Any: ...

class ReferenceType(Generic[_T]):
    __callback__: Callable[[ReferenceType[_T]], Any]
    def __init__(self, o: _T, callback: Callable[[ReferenceType[_T]], Any] | None = ...) -> None: ...
    def __call__(self) -> _T | None: ...
    def __hash__(self) -> int: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

ref = ReferenceType

def getweakrefcount(__object: Any) -> int: ...
def getweakrefs(object: Any) -> list[Any]: ...
@overload
def proxy(object: _C, callback: Callable[[_C], Any] | None = ...) -> CallableProxyType[_C]: ...

# Return CallableProxyType if object is callable, ProxyType otherwise
@overload
def proxy(object: _T, callback: Callable[[_T], Any] | None = ...) -> Any: ...
