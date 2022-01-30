from typing import Dict, List, Optional, Text, Tuple

class NetrcParseError(Exception):
    filename: str | None
    lineno: int | None
    msg: str
    def __init__(self, msg: str, filename: Text | None = ..., lineno: int | None = ...) -> None: ...

# (login, account, password) tuple
_NetrcTuple = Tuple[str, Optional[str], Optional[str]]

class netrc:
    hosts: Dict[str, _NetrcTuple]
    macros: Dict[str, List[str]]
    def __init__(self, file: Text | None = ...) -> None: ...
    def authenticators(self, host: str) -> _NetrcTuple | None: ...
