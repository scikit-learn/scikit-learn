from typing import Any, Callable

class JSONDecodeError(ValueError):
    msg: str
    doc: str
    pos: int
    lineno: int
    colno: int
    def __init__(self, msg: str, doc: str, pos: int) -> None: ...

class JSONDecoder:
    object_hook: Callable[[dict[str, Any]], Any]
    parse_float: Callable[[str], Any]
    parse_int: Callable[[str], Any]
    parse_constant: Callable[[str], Any]
    strict: bool
    object_pairs_hook: Callable[[list[tuple[str, Any]]], Any]
    def __init__(
        self,
        *,
        object_hook: Callable[[dict[str, Any]], Any] | None = ...,
        parse_float: Callable[[str], Any] | None = ...,
        parse_int: Callable[[str], Any] | None = ...,
        parse_constant: Callable[[str], Any] | None = ...,
        strict: bool = ...,
        object_pairs_hook: Callable[[list[tuple[str, Any]]], Any] | None = ...,
    ) -> None: ...
    def decode(self, s: str, _w: Callable[..., Any] = ...) -> Any: ...  # _w is undocumented
    def raw_decode(self, s: str, idx: int = ...) -> tuple[Any, int]: ...
