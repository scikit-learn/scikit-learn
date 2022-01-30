from typing import Any, Callable, NamedTuple, Tuple

class Event(NamedTuple):
    time: float
    priority: Any
    action: Callable[..., Any]
    argument: Tuple[Any, ...]
    kwargs: dict[str, Any]

class scheduler:
    def __init__(self, timefunc: Callable[[], float] = ..., delayfunc: Callable[[float], None] = ...) -> None: ...
    def enterabs(
        self,
        time: float,
        priority: Any,
        action: Callable[..., Any],
        argument: Tuple[Any, ...] = ...,
        kwargs: dict[str, Any] = ...,
    ) -> Event: ...
    def enter(
        self,
        delay: float,
        priority: Any,
        action: Callable[..., Any],
        argument: Tuple[Any, ...] = ...,
        kwargs: dict[str, Any] = ...,
    ) -> Event: ...
    def run(self, blocking: bool = ...) -> float | None: ...
    def cancel(self, event: Event) -> None: ...
    def empty(self) -> bool: ...
    @property
    def queue(self) -> list[Event]: ...
