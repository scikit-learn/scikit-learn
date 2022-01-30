import sys
from typing import Any, AnyStr, Callable, Generic, Iterable, Iterator, NamedTuple, Sequence, TypeVar, Union, overload

if sys.version_info >= (3, 9):
    from types import GenericAlias

_T = TypeVar("_T")
_JunkCallback = Union[Callable[[str], bool], Callable[[str], bool]]

class Match(NamedTuple):
    a: int
    b: int
    size: int

class SequenceMatcher(Generic[_T]):
    def __init__(
        self, isjunk: Callable[[_T], bool] | None = ..., a: Sequence[_T] = ..., b: Sequence[_T] = ..., autojunk: bool = ...
    ) -> None: ...
    def set_seqs(self, a: Sequence[_T], b: Sequence[_T]) -> None: ...
    def set_seq1(self, a: Sequence[_T]) -> None: ...
    def set_seq2(self, b: Sequence[_T]) -> None: ...
    if sys.version_info >= (3, 9):
        def find_longest_match(self, alo: int = ..., ahi: int | None = ..., blo: int = ..., bhi: int | None = ...) -> Match: ...
    else:
        def find_longest_match(self, alo: int, ahi: int, blo: int, bhi: int) -> Match: ...
    def get_matching_blocks(self) -> list[Match]: ...
    def get_opcodes(self) -> list[tuple[str, int, int, int, int]]: ...
    def get_grouped_opcodes(self, n: int = ...) -> Iterable[list[tuple[str, int, int, int, int]]]: ...
    def ratio(self) -> float: ...
    def quick_ratio(self) -> float: ...
    def real_quick_ratio(self) -> float: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

# mypy thinks the signatures of the overloads overlap, but the types still work fine
@overload
def get_close_matches(  # type: ignore
    word: AnyStr, possibilities: Iterable[AnyStr], n: int = ..., cutoff: float = ...
) -> list[AnyStr]: ...
@overload
def get_close_matches(
    word: Sequence[_T], possibilities: Iterable[Sequence[_T]], n: int = ..., cutoff: float = ...
) -> list[Sequence[_T]]: ...

class Differ:
    def __init__(self, linejunk: _JunkCallback | None = ..., charjunk: _JunkCallback | None = ...) -> None: ...
    def compare(self, a: Sequence[str], b: Sequence[str]) -> Iterator[str]: ...

def IS_LINE_JUNK(line: str, pat: Any = ...) -> bool: ...  # pat is undocumented
def IS_CHARACTER_JUNK(ch: str, ws: str = ...) -> bool: ...  # ws is undocumented
def unified_diff(
    a: Sequence[str],
    b: Sequence[str],
    fromfile: str = ...,
    tofile: str = ...,
    fromfiledate: str = ...,
    tofiledate: str = ...,
    n: int = ...,
    lineterm: str = ...,
) -> Iterator[str]: ...
def context_diff(
    a: Sequence[str],
    b: Sequence[str],
    fromfile: str = ...,
    tofile: str = ...,
    fromfiledate: str = ...,
    tofiledate: str = ...,
    n: int = ...,
    lineterm: str = ...,
) -> Iterator[str]: ...
def ndiff(
    a: Sequence[str], b: Sequence[str], linejunk: _JunkCallback | None = ..., charjunk: _JunkCallback | None = ...
) -> Iterator[str]: ...

class HtmlDiff(object):
    def __init__(
        self,
        tabsize: int = ...,
        wrapcolumn: int | None = ...,
        linejunk: _JunkCallback | None = ...,
        charjunk: _JunkCallback | None = ...,
    ) -> None: ...
    def make_file(
        self,
        fromlines: Sequence[str],
        tolines: Sequence[str],
        fromdesc: str = ...,
        todesc: str = ...,
        context: bool = ...,
        numlines: int = ...,
        *,
        charset: str = ...,
    ) -> str: ...
    def make_table(
        self,
        fromlines: Sequence[str],
        tolines: Sequence[str],
        fromdesc: str = ...,
        todesc: str = ...,
        context: bool = ...,
        numlines: int = ...,
    ) -> str: ...

def restore(delta: Iterable[str], which: int) -> Iterator[str]: ...
def diff_bytes(
    dfunc: Callable[[Sequence[str], Sequence[str], str, str, str, str, int, str], Iterator[str]],
    a: Sequence[bytes],
    b: Sequence[bytes],
    fromfile: bytes = ...,
    tofile: bytes = ...,
    fromfiledate: bytes = ...,
    tofiledate: bytes = ...,
    n: int = ...,
    lineterm: bytes = ...,
) -> Iterator[bytes]: ...
