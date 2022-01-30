import sys
from _typeshed import SupportsWrite
from types import FrameType, TracebackType
from typing import IO, Any, Generator, Iterable, Iterator, List, Mapping, Optional, Tuple, Type

_PT = Tuple[str, int, str, Optional[str]]

def print_tb(tb: TracebackType | None, limit: int | None = ..., file: IO[str] | None = ...) -> None: ...

if sys.version_info >= (3, 10):
    def print_exception(
        __exc: Type[BaseException] | None,
        value: BaseException | None = ...,
        tb: TracebackType | None = ...,
        limit: int | None = ...,
        file: IO[str] | None = ...,
        chain: bool = ...,
    ) -> None: ...

else:
    def print_exception(
        etype: Type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
        limit: int | None = ...,
        file: IO[str] | None = ...,
        chain: bool = ...,
    ) -> None: ...

def print_exc(limit: int | None = ..., file: IO[str] | None = ..., chain: bool = ...) -> None: ...
def print_last(limit: int | None = ..., file: IO[str] | None = ..., chain: bool = ...) -> None: ...
def print_stack(f: FrameType | None = ..., limit: int | None = ..., file: IO[str] | None = ...) -> None: ...
def extract_tb(tb: TracebackType | None, limit: int | None = ...) -> StackSummary: ...
def extract_stack(f: FrameType | None = ..., limit: int | None = ...) -> StackSummary: ...
def format_list(extracted_list: list[FrameSummary]) -> list[str]: ...

# undocumented
def print_list(extracted_list: list[FrameSummary], file: SupportsWrite[str] | None = ...) -> None: ...

if sys.version_info >= (3, 10):
    def format_exception_only(__exc: Type[BaseException] | None, value: BaseException | None = ...) -> list[str]: ...

else:
    def format_exception_only(etype: Type[BaseException] | None, value: BaseException | None) -> list[str]: ...

if sys.version_info >= (3, 10):
    def format_exception(
        __exc: Type[BaseException] | None,
        value: BaseException | None = ...,
        tb: TracebackType | None = ...,
        limit: int | None = ...,
        chain: bool = ...,
    ) -> list[str]: ...

else:
    def format_exception(
        etype: Type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
        limit: int | None = ...,
        chain: bool = ...,
    ) -> list[str]: ...

def format_exc(limit: int | None = ..., chain: bool = ...) -> str: ...
def format_tb(tb: TracebackType | None, limit: int | None = ...) -> list[str]: ...
def format_stack(f: FrameType | None = ..., limit: int | None = ...) -> list[str]: ...
def clear_frames(tb: TracebackType) -> None: ...
def walk_stack(f: FrameType | None) -> Iterator[tuple[FrameType, int]]: ...
def walk_tb(tb: TracebackType | None) -> Iterator[tuple[FrameType, int]]: ...

class TracebackException:
    __cause__: TracebackException
    __context__: TracebackException
    __suppress_context__: bool
    stack: StackSummary
    exc_type: Type[BaseException]
    filename: str
    lineno: int
    text: str
    offset: int
    msg: str
    if sys.version_info >= (3, 10):
        def __init__(
            self,
            exc_type: Type[BaseException],
            exc_value: BaseException,
            exc_traceback: TracebackType,
            *,
            limit: int | None = ...,
            lookup_lines: bool = ...,
            capture_locals: bool = ...,
            compact: bool = ...,
            _seen: set[int] | None = ...,
        ) -> None: ...
        @classmethod
        def from_exception(
            cls,
            exc: BaseException,
            *,
            limit: int | None = ...,
            lookup_lines: bool = ...,
            capture_locals: bool = ...,
            compact: bool = ...,
        ) -> TracebackException: ...
    else:
        def __init__(
            self,
            exc_type: Type[BaseException],
            exc_value: BaseException,
            exc_traceback: TracebackType,
            *,
            limit: int | None = ...,
            lookup_lines: bool = ...,
            capture_locals: bool = ...,
            _seen: set[int] | None = ...,
        ) -> None: ...
        @classmethod
        def from_exception(
            cls, exc: BaseException, *, limit: int | None = ..., lookup_lines: bool = ..., capture_locals: bool = ...
        ) -> TracebackException: ...
    def format(self, *, chain: bool = ...) -> Generator[str, None, None]: ...
    def format_exception_only(self) -> Generator[str, None, None]: ...

class FrameSummary(Iterable[Any]):
    filename: str
    lineno: int
    name: str
    line: str
    locals: dict[str, str] | None
    def __init__(
        self,
        filename: str,
        lineno: int,
        name: str,
        *,
        lookup_line: bool = ...,
        locals: Mapping[str, str] | None = ...,
        line: str | None = ...,
    ) -> None: ...
    # TODO: more precise typing for __getitem__ and __iter__,
    # for a namedtuple-like view on (filename, lineno, name, str).
    def __getitem__(self, i: int) -> Any: ...
    def __iter__(self) -> Iterator[Any]: ...

class StackSummary(List[FrameSummary]):
    @classmethod
    def extract(
        cls,
        frame_gen: Iterable[tuple[FrameType, int]],
        *,
        limit: int | None = ...,
        lookup_lines: bool = ...,
        capture_locals: bool = ...,
    ) -> StackSummary: ...
    @classmethod
    def from_list(cls, a_list: list[_PT]) -> StackSummary: ...
    def format(self) -> list[str]: ...
