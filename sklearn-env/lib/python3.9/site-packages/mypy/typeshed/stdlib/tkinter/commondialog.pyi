from typing import Any, ClassVar, Mapping

class Dialog:
    command: ClassVar[str | None]
    master: Any | None
    options: Mapping[str, Any]
    def __init__(self, master: Any | None = ..., **options) -> None: ...
    def show(self, **options) -> Any: ...
