from typing import ClassVar, Literal

from .. import fixer_base

class FixSysExc(fixer_base.BaseFix):
    exc_info: ClassVar[list[str]]
    BM_compatible: ClassVar[Literal[True]]
    PATTERN: ClassVar[str]
    def transform(self, node, results): ...
