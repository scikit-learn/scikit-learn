from typing import ClassVar, Literal

from _typeshed import Incomplete

from .. import fixer_base

iter_exempt: set[str]

class FixDict(fixer_base.BaseFix):
    BM_compatible: ClassVar[Literal[True]]
    PATTERN: ClassVar[str]
    def transform(self, node, results): ...
    P1: ClassVar[str]
    p1: ClassVar[Incomplete]
    P2: ClassVar[str]
    p2: ClassVar[Incomplete]
    def in_special_context(self, node, isiter): ...
