from typing import Any, Callable, Hashable, Optional, SupportsInt, Tuple, TypeVar, Union

_TypeT = TypeVar("_TypeT", bound=type)
_Reduce = Union[Tuple[Callable[..., _TypeT], Tuple[Any, ...]], Tuple[Callable[..., _TypeT], Tuple[Any, ...], Optional[Any]]]

__all__: list[str]

def pickle(
    ob_type: _TypeT,
    pickle_function: Callable[[_TypeT], str | _Reduce[_TypeT]],
    constructor_ob: Callable[[_Reduce[_TypeT]], _TypeT] | None = ...,
) -> None: ...
def constructor(object: Callable[[_Reduce[_TypeT]], _TypeT]) -> None: ...
def add_extension(module: Hashable, name: Hashable, code: SupportsInt) -> None: ...
def remove_extension(module: Hashable, name: Hashable, code: int) -> None: ...
def clear_extension_cache() -> None: ...

dispatch_table: dict[type, Callable[[type], str | _Reduce[type]]]  # undocumented
