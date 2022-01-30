from typing import Any, Dict, Union

_Text = Union[str, unicode]

class Completer:
    def __init__(self, namespace: Dict[str, Any] | None = ...) -> None: ...
    def complete(self, text: _Text, state: int) -> str | None: ...
