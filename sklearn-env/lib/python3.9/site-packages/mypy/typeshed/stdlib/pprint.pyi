import sys
from typing import IO, Any

if sys.version_info >= (3, 10):
    def pformat(
        object: object,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
        underscore_numbers: bool = ...,
    ) -> str: ...

elif sys.version_info >= (3, 8):
    def pformat(
        object: object,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
    ) -> str: ...

else:
    def pformat(object: object, indent: int = ..., width: int = ..., depth: int | None = ..., *, compact: bool = ...) -> str: ...

if sys.version_info >= (3, 10):
    def pp(
        object: object,
        stream: IO[str] | None = ...,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
        underscore_numbers: bool = ...,
    ) -> None: ...

elif sys.version_info >= (3, 8):
    def pp(
        object: object,
        stream: IO[str] | None = ...,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
    ) -> None: ...

if sys.version_info >= (3, 10):
    def pprint(
        object: object,
        stream: IO[str] | None = ...,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
        underscore_numbers: bool = ...,
    ) -> None: ...

elif sys.version_info >= (3, 8):
    def pprint(
        object: object,
        stream: IO[str] | None = ...,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
        sort_dicts: bool = ...,
    ) -> None: ...

else:
    def pprint(
        object: object,
        stream: IO[str] | None = ...,
        indent: int = ...,
        width: int = ...,
        depth: int | None = ...,
        *,
        compact: bool = ...,
    ) -> None: ...

def isreadable(object: object) -> bool: ...
def isrecursive(object: object) -> bool: ...
def saferepr(object: object) -> str: ...

class PrettyPrinter:
    if sys.version_info >= (3, 10):
        def __init__(
            self,
            indent: int = ...,
            width: int = ...,
            depth: int | None = ...,
            stream: IO[str] | None = ...,
            *,
            compact: bool = ...,
            sort_dicts: bool = ...,
            underscore_numbers: bool = ...,
        ) -> None: ...
    elif sys.version_info >= (3, 8):
        def __init__(
            self,
            indent: int = ...,
            width: int = ...,
            depth: int | None = ...,
            stream: IO[str] | None = ...,
            *,
            compact: bool = ...,
            sort_dicts: bool = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            indent: int = ...,
            width: int = ...,
            depth: int | None = ...,
            stream: IO[str] | None = ...,
            *,
            compact: bool = ...,
        ) -> None: ...
    def pformat(self, object: object) -> str: ...
    def pprint(self, object: object) -> None: ...
    def isreadable(self, object: object) -> bool: ...
    def isrecursive(self, object: object) -> bool: ...
    def format(self, object: object, context: dict[int, Any], maxlevels: int, level: int) -> tuple[str, bool, bool]: ...
