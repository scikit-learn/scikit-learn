import sys
import types
from collections.abc import Callable, Coroutine
from typing import Any, TypeVar
from typing_extensions import TypeGuard

_F = TypeVar("_F", bound=Callable[..., Any])

def coroutine(func: _F) -> _F: ...
def iscoroutinefunction(func: object) -> bool: ...

if sys.version_info < (3, 8):
    def iscoroutine(obj: object) -> TypeGuard[types.GeneratorType[Any, Any, Any] | Coroutine[Any, Any, Any]]: ...

else:
    def iscoroutine(obj: object) -> TypeGuard[Coroutine[Any, Any, Any]]: ...
