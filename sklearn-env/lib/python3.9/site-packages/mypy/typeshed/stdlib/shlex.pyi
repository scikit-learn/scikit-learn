import sys
from typing import Any, Iterable, TextIO, TypeVar

def split(s: str, comments: bool = ..., posix: bool = ...) -> list[str]: ...

if sys.version_info >= (3, 8):
    def join(split_command: Iterable[str]) -> str: ...

def quote(s: str) -> str: ...

_SLT = TypeVar("_SLT", bound=shlex)

class shlex(Iterable[str]):
    commenters: str
    wordchars: str
    whitespace: str
    escape: str
    quotes: str
    escapedquotes: str
    whitespace_split: bool
    infile: str
    instream: TextIO
    source: str
    debug: int
    lineno: int
    token: str
    eof: str
    punctuation_chars: str
    def __init__(
        self,
        instream: str | TextIO | None = ...,
        infile: str | None = ...,
        posix: bool = ...,
        punctuation_chars: bool | str = ...,
    ) -> None: ...
    def get_token(self) -> str: ...
    def push_token(self, tok: str) -> None: ...
    def read_token(self) -> str: ...
    def sourcehook(self, newfile: str) -> tuple[str, TextIO]: ...
    # TODO argument types
    def push_source(self, newstream: Any, newfile: Any = ...) -> None: ...
    def pop_source(self) -> None: ...
    def error_leader(self, infile: str | None = ..., lineno: int | None = ...) -> None: ...
    def __iter__(self: _SLT) -> _SLT: ...
    def __next__(self) -> str: ...
