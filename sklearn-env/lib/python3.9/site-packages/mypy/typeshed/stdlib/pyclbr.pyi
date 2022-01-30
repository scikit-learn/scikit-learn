import sys
from typing import Sequence

class Class:
    module: str
    name: str
    super: list[Class | str] | None
    methods: dict[str, int]
    file: int
    lineno: int

    if sys.version_info >= (3, 7):
        def __init__(
            self, module: str, name: str, super: list[Class | str] | None, file: str, lineno: int, parent: Class | None = ...
        ) -> None: ...
    else:
        def __init__(self, module: str, name: str, super: list[Class | str] | None, file: str, lineno: int) -> None: ...

class Function:
    module: str
    name: str
    file: int
    lineno: int

    if sys.version_info >= (3, 7):
        def __init__(self, module: str, name: str, file: str, lineno: int, parent: Function | None = ...) -> None: ...
    else:
        def __init__(self, module: str, name: str, file: str, lineno: int) -> None: ...

def readmodule(module: str, path: Sequence[str] | None = ...) -> dict[str, Class]: ...
def readmodule_ex(module: str, path: Sequence[str] | None = ...) -> dict[str, Class | Function | list[str]]: ...
