from typing import ClassVar, Literal

from .. import fixer_base

class FixItertools(fixer_base.BaseFix):
    BM_compatible: ClassVar[Literal[True]]
    it_funcs: str
    PATTERN: ClassVar[str]
    def transform(self, node, results) -> None: ...
