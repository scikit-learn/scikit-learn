import sys
from _typeshed import Self
from typing import Any, Callable, ContextManager, Generic, Iterable, Iterator, List, Mapping, TypeVar

if sys.version_info >= (3, 9):
    from types import GenericAlias

_PT = TypeVar("_PT", bound=Pool)
_S = TypeVar("_S")
_T = TypeVar("_T")

class ApplyResult(Generic[_T]):
    if sys.version_info >= (3, 8):
        def __init__(
            self, pool: Pool, callback: Callable[[_T], None] | None, error_callback: Callable[[BaseException], None] | None
        ) -> None: ...
    else:
        def __init__(
            self,
            cache: dict[int, ApplyResult[Any]],
            callback: Callable[[_T], None] | None,
            error_callback: Callable[[BaseException], None] | None,
        ) -> None: ...
    def get(self, timeout: float | None = ...) -> _T: ...
    def wait(self, timeout: float | None = ...) -> None: ...
    def ready(self) -> bool: ...
    def successful(self) -> bool: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

# alias created during issue #17805
AsyncResult = ApplyResult

class MapResult(ApplyResult[List[_T]]):
    if sys.version_info >= (3, 8):
        def __init__(
            self,
            pool: Pool,
            chunksize: int,
            length: int,
            callback: Callable[[list[_T]], None] | None,
            error_callback: Callable[[BaseException], None] | None,
        ) -> None: ...
    else:
        def __init__(
            self,
            cache: dict[int, ApplyResult[Any]],
            chunksize: int,
            length: int,
            callback: Callable[[list[_T]], None] | None,
            error_callback: Callable[[BaseException], None] | None,
        ) -> None: ...

class IMapIterator(Iterator[_T]):
    if sys.version_info >= (3, 8):
        def __init__(self, pool: Pool) -> None: ...
    else:
        def __init__(self, cache: dict[int, IMapIterator[Any]]) -> None: ...
    def __iter__(self: _S) -> _S: ...
    def next(self, timeout: float | None = ...) -> _T: ...
    def __next__(self, timeout: float | None = ...) -> _T: ...

class IMapUnorderedIterator(IMapIterator[_T]): ...

class Pool(ContextManager[Pool]):
    def __init__(
        self,
        processes: int | None = ...,
        initializer: Callable[..., None] | None = ...,
        initargs: Iterable[Any] = ...,
        maxtasksperchild: int | None = ...,
        context: Any | None = ...,
    ) -> None: ...
    def apply(self, func: Callable[..., _T], args: Iterable[Any] = ..., kwds: Mapping[str, Any] = ...) -> _T: ...
    def apply_async(
        self,
        func: Callable[..., _T],
        args: Iterable[Any] = ...,
        kwds: Mapping[str, Any] = ...,
        callback: Callable[[_T], None] | None = ...,
        error_callback: Callable[[BaseException], None] | None = ...,
    ) -> AsyncResult[_T]: ...
    def map(self, func: Callable[[_S], _T], iterable: Iterable[_S], chunksize: int | None = ...) -> list[_T]: ...
    def map_async(
        self,
        func: Callable[[_S], _T],
        iterable: Iterable[_S],
        chunksize: int | None = ...,
        callback: Callable[[_T], None] | None = ...,
        error_callback: Callable[[BaseException], None] | None = ...,
    ) -> MapResult[_T]: ...
    def imap(self, func: Callable[[_S], _T], iterable: Iterable[_S], chunksize: int | None = ...) -> IMapIterator[_T]: ...
    def imap_unordered(
        self, func: Callable[[_S], _T], iterable: Iterable[_S], chunksize: int | None = ...
    ) -> IMapIterator[_T]: ...
    def starmap(self, func: Callable[..., _T], iterable: Iterable[Iterable[Any]], chunksize: int | None = ...) -> list[_T]: ...
    def starmap_async(
        self,
        func: Callable[..., _T],
        iterable: Iterable[Iterable[Any]],
        chunksize: int | None = ...,
        callback: Callable[[_T], None] | None = ...,
        error_callback: Callable[[BaseException], None] | None = ...,
    ) -> AsyncResult[list[_T]]: ...
    def close(self) -> None: ...
    def terminate(self) -> None: ...
    def join(self) -> None: ...
    def __enter__(self: Self) -> Self: ...

class ThreadPool(Pool, ContextManager[ThreadPool]):
    def __init__(
        self, processes: int | None = ..., initializer: Callable[..., Any] | None = ..., initargs: Iterable[Any] = ...
    ) -> None: ...

# undocumented
if sys.version_info >= (3, 8):
    INIT: str
    RUN: str
    CLOSE: str
    TERMINATE: str
else:
    RUN: int
    CLOSE: int
    TERMINATE: int
