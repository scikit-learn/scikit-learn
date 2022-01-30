import concurrent.futures
import sys
from collections.abc import Awaitable, Generator, Iterable, Iterator
from types import FrameType
from typing import Any, Generic, Optional, TextIO, TypeVar, Union, overload
from typing_extensions import Literal

from .events import AbstractEventLoop
from .futures import Future

if sys.version_info >= (3, 9):
    from types import GenericAlias

_T = TypeVar("_T")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_T3 = TypeVar("_T3")
_T4 = TypeVar("_T4")
_T5 = TypeVar("_T5")
_FT = TypeVar("_FT", bound=Future[Any])
_FutureT = Union[Future[_T], Generator[Any, None, _T], Awaitable[_T]]
_TaskYieldType = Optional[Future[object]]

FIRST_EXCEPTION: str
FIRST_COMPLETED: str
ALL_COMPLETED: str

if sys.version_info >= (3, 10):
    def as_completed(fs: Iterable[_FutureT[_T]], *, timeout: float | None = ...) -> Iterator[Future[_T]]: ...

else:
    def as_completed(
        fs: Iterable[_FutureT[_T]], *, loop: AbstractEventLoop | None = ..., timeout: float | None = ...
    ) -> Iterator[Future[_T]]: ...

@overload
def ensure_future(coro_or_future: _FT, *, loop: AbstractEventLoop | None = ...) -> _FT: ...  # type: ignore
@overload
def ensure_future(coro_or_future: Awaitable[_T], *, loop: AbstractEventLoop | None = ...) -> Task[_T]: ...

# Prior to Python 3.7 'async' was an alias for 'ensure_future'.
# It became a keyword in 3.7.

# `gather()` actually returns a list with length equal to the number
# of tasks passed; however, Tuple is used similar to the annotation for
# zip() because typing does not support variadic type variables.  See
# typing PR #1550 for discussion.
if sys.version_info >= (3, 10):
    @overload
    def gather(coro_or_future1: _FutureT[_T1], *, return_exceptions: Literal[False] = ...) -> Future[tuple[_T1]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1], coro_or_future2: _FutureT[_T2], *, return_exceptions: Literal[False] = ...
    ) -> Future[tuple[_T1, _T2]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        *,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        *,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3, _T4]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        coro_or_future5: _FutureT[_T5],
        *,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3, _T4, _T5]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[Any],
        coro_or_future2: _FutureT[Any],
        coro_or_future3: _FutureT[Any],
        coro_or_future4: _FutureT[Any],
        coro_or_future5: _FutureT[Any],
        coro_or_future6: _FutureT[Any],
        *coros_or_futures: _FutureT[Any],
        return_exceptions: bool = ...,
    ) -> Future[list[Any]]: ...
    @overload
    def gather(coro_or_future1: _FutureT[_T1], *, return_exceptions: bool = ...) -> Future[tuple[_T1 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1], coro_or_future2: _FutureT[_T2], *, return_exceptions: bool = ...
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        *,
        return_exceptions: bool = ...,
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        *,
        return_exceptions: bool = ...,
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException, _T4 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        coro_or_future5: _FutureT[_T5],
        *,
        return_exceptions: bool = ...,
    ) -> Future[
        tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException, _T4 | BaseException, _T5 | BaseException]
    ]: ...

else:
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1], *, loop: AbstractEventLoop | None = ..., return_exceptions: Literal[False] = ...
    ) -> Future[tuple[_T1]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3, _T4]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        coro_or_future5: _FutureT[_T5],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: Literal[False] = ...,
    ) -> Future[tuple[_T1, _T2, _T3, _T4, _T5]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[Any],
        coro_or_future2: _FutureT[Any],
        coro_or_future3: _FutureT[Any],
        coro_or_future4: _FutureT[Any],
        coro_or_future5: _FutureT[Any],
        coro_or_future6: _FutureT[Any],
        *coros_or_futures: _FutureT[Any],
        loop: AbstractEventLoop | None = ...,
        return_exceptions: bool = ...,
    ) -> Future[list[Any]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1], *, loop: AbstractEventLoop | None = ..., return_exceptions: bool = ...
    ) -> Future[tuple[_T1 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: bool = ...,
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: bool = ...,
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: bool = ...,
    ) -> Future[tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException, _T4 | BaseException]]: ...
    @overload
    def gather(
        coro_or_future1: _FutureT[_T1],
        coro_or_future2: _FutureT[_T2],
        coro_or_future3: _FutureT[_T3],
        coro_or_future4: _FutureT[_T4],
        coro_or_future5: _FutureT[_T5],
        *,
        loop: AbstractEventLoop | None = ...,
        return_exceptions: bool = ...,
    ) -> Future[
        tuple[_T1 | BaseException, _T2 | BaseException, _T3 | BaseException, _T4 | BaseException, _T5 | BaseException]
    ]: ...

def run_coroutine_threadsafe(coro: _FutureT[_T], loop: AbstractEventLoop) -> concurrent.futures.Future[_T]: ...

if sys.version_info >= (3, 10):
    def shield(arg: _FutureT[_T]) -> Future[_T]: ...
    def sleep(delay: float, result: _T = ...) -> Future[_T]: ...
    @overload
    def wait(fs: Iterable[_FT], *, timeout: float | None = ..., return_when: str = ...) -> Future[tuple[set[_FT], set[_FT]]]: ...  # type: ignore
    @overload
    def wait(
        fs: Iterable[Awaitable[_T]], *, timeout: float | None = ..., return_when: str = ...
    ) -> Future[tuple[set[Task[_T]], set[Task[_T]]]]: ...
    def wait_for(fut: _FutureT[_T], timeout: float | None) -> Future[_T]: ...

else:
    def shield(arg: _FutureT[_T], *, loop: AbstractEventLoop | None = ...) -> Future[_T]: ...
    def sleep(delay: float, result: _T = ..., *, loop: AbstractEventLoop | None = ...) -> Future[_T]: ...
    @overload
    def wait(fs: Iterable[_FT], *, loop: AbstractEventLoop | None = ..., timeout: float | None = ..., return_when: str = ...) -> Future[tuple[set[_FT], set[_FT]]]: ...  # type: ignore
    @overload
    def wait(
        fs: Iterable[Awaitable[_T]], *, loop: AbstractEventLoop | None = ..., timeout: float | None = ..., return_when: str = ...
    ) -> Future[tuple[set[Task[_T]], set[Task[_T]]]]: ...
    def wait_for(fut: _FutureT[_T], timeout: float | None, *, loop: AbstractEventLoop | None = ...) -> Future[_T]: ...

class Task(Future[_T], Generic[_T]):
    if sys.version_info >= (3, 8):
        def __init__(
            self,
            coro: Generator[_TaskYieldType, None, _T] | Awaitable[_T],
            *,
            loop: AbstractEventLoop = ...,
            name: str | None = ...,
        ) -> None: ...
    else:
        def __init__(
            self, coro: Generator[_TaskYieldType, None, _T] | Awaitable[_T], *, loop: AbstractEventLoop = ...
        ) -> None: ...
    def __repr__(self) -> str: ...
    if sys.version_info >= (3, 8):
        def get_coro(self) -> Generator[_TaskYieldType, None, _T] | Awaitable[_T]: ...
        def get_name(self) -> str: ...
        def set_name(self, __value: object) -> None: ...
    def get_stack(self, *, limit: int | None = ...) -> list[FrameType]: ...
    def print_stack(self, *, limit: int | None = ..., file: TextIO | None = ...) -> None: ...
    if sys.version_info >= (3, 9):
        def cancel(self, msg: Any | None = ...) -> bool: ...
    else:
        def cancel(self) -> bool: ...
    if sys.version_info < (3, 9):
        @classmethod
        def current_task(cls, loop: AbstractEventLoop | None = ...) -> Task[Any] | None: ...
        @classmethod
        def all_tasks(cls, loop: AbstractEventLoop | None = ...) -> set[Task[Any]]: ...
    if sys.version_info < (3, 7):
        def _wakeup(self, fut: Future[Any]) -> None: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

if sys.version_info >= (3, 7):
    def all_tasks(loop: AbstractEventLoop | None = ...) -> set[Task[Any]]: ...
    if sys.version_info >= (3, 8):
        def create_task(coro: Generator[_TaskYieldType, None, _T] | Awaitable[_T], *, name: str | None = ...) -> Task[_T]: ...
    else:
        def create_task(coro: Generator[_TaskYieldType, None, _T] | Awaitable[_T]) -> Task[_T]: ...
    def current_task(loop: AbstractEventLoop | None = ...) -> Task[Any] | None: ...
