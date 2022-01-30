from types import TracebackType
from typing import IO, Any, Callable, ContextManager, Iterable, Iterator, Optional, Type, TypeVar
from typing_extensions import Protocol

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)
_T_io = TypeVar("_T_io", bound=Optional[IO[str]])
_F = TypeVar("_F", bound=Callable[..., Any])

_ExitFunc = Callable[[Optional[Type[BaseException]], Optional[BaseException], Optional[TracebackType]], bool]
_CM_EF = TypeVar("_CM_EF", ContextManager[Any], _ExitFunc)

class GeneratorContextManager(ContextManager[_T_co]):
    def __call__(self, func: _F) -> _F: ...

def contextmanager(func: Callable[..., Iterator[_T]]) -> Callable[..., ContextManager[_T]]: ...
def nested(*mgr: ContextManager[Any]) -> ContextManager[Iterable[Any]]: ...

class _SupportsClose(Protocol):
    def close(self) -> None: ...

_SupportsCloseT = TypeVar("_SupportsCloseT", bound=_SupportsClose)

class closing(ContextManager[_SupportsCloseT]):
    def __init__(self, thing: _SupportsCloseT) -> None: ...
