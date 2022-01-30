from _typeshed import StrPath
from typing import Dict, List, Optional, Tuple, TypeVar

_P = TypeVar("_P")
_Label = Tuple[int, Optional[str]]
_DFA = List[List[Tuple[int, int]]]
_DFAS = Tuple[_DFA, Dict[int, int]]

class Grammar:
    symbol2number: dict[str, int]
    number2symbol: dict[int, str]
    states: list[_DFA]
    dfas: dict[int, _DFAS]
    labels: list[_Label]
    keywords: dict[str, int]
    tokens: dict[int, int]
    symbol2label: dict[str, int]
    start: int
    def __init__(self) -> None: ...
    def dump(self, filename: StrPath) -> None: ...
    def load(self, filename: StrPath) -> None: ...
    def copy(self: _P) -> _P: ...
    def report(self) -> None: ...

opmap_raw: str
opmap: dict[str, str]
