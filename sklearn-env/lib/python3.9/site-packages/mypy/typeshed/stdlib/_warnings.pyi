from typing import Any, Type, overload

_defaultaction: str
_onceregistry: dict[Any, Any]
filters: list[tuple[str, str | None, Type[Warning], str | None, int]]

@overload
def warn(message: str, category: Type[Warning] | None = ..., stacklevel: int = ..., source: Any | None = ...) -> None: ...
@overload
def warn(message: Warning, category: Any = ..., stacklevel: int = ..., source: Any | None = ...) -> None: ...
@overload
def warn_explicit(
    message: str,
    category: Type[Warning],
    filename: str,
    lineno: int,
    module: str | None = ...,
    registry: dict[str | tuple[str, Type[Warning], int], int] | None = ...,
    module_globals: dict[str, Any] | None = ...,
    source: Any | None = ...,
) -> None: ...
@overload
def warn_explicit(
    message: Warning,
    category: Any,
    filename: str,
    lineno: int,
    module: str | None = ...,
    registry: dict[str | tuple[str, Type[Warning], int], int] | None = ...,
    module_globals: dict[str, Any] | None = ...,
    source: Any | None = ...,
) -> None: ...
