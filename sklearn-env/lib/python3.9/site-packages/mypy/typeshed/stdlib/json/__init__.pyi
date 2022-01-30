from _typeshed import SupportsRead
from typing import IO, Any, Callable, Type

from .decoder import JSONDecodeError as JSONDecodeError, JSONDecoder as JSONDecoder
from .encoder import JSONEncoder as JSONEncoder

def dumps(
    obj: Any,
    *,
    skipkeys: bool = ...,
    ensure_ascii: bool = ...,
    check_circular: bool = ...,
    allow_nan: bool = ...,
    cls: Type[JSONEncoder] | None = ...,
    indent: None | int | str = ...,
    separators: tuple[str, str] | None = ...,
    default: Callable[[Any], Any] | None = ...,
    sort_keys: bool = ...,
    **kwds: Any,
) -> str: ...
def dump(
    obj: Any,
    fp: IO[str],
    *,
    skipkeys: bool = ...,
    ensure_ascii: bool = ...,
    check_circular: bool = ...,
    allow_nan: bool = ...,
    cls: Type[JSONEncoder] | None = ...,
    indent: None | int | str = ...,
    separators: tuple[str, str] | None = ...,
    default: Callable[[Any], Any] | None = ...,
    sort_keys: bool = ...,
    **kwds: Any,
) -> None: ...
def loads(
    s: str | bytes,
    *,
    cls: Type[JSONDecoder] | None = ...,
    object_hook: Callable[[dict[Any, Any]], Any] | None = ...,
    parse_float: Callable[[str], Any] | None = ...,
    parse_int: Callable[[str], Any] | None = ...,
    parse_constant: Callable[[str], Any] | None = ...,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = ...,
    **kwds: Any,
) -> Any: ...
def load(
    fp: SupportsRead[str | bytes],
    *,
    cls: Type[JSONDecoder] | None = ...,
    object_hook: Callable[[dict[Any, Any]], Any] | None = ...,
    parse_float: Callable[[str], Any] | None = ...,
    parse_int: Callable[[str], Any] | None = ...,
    parse_constant: Callable[[str], Any] | None = ...,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = ...,
    **kwds: Any,
) -> Any: ...
def detect_encoding(b: bytes) -> str: ...  # undocumented
