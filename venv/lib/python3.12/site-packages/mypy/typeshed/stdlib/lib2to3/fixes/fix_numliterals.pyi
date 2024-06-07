from typing import ClassVar, Literal

from .. import fixer_base

class FixNumliterals(fixer_base.BaseFix):
    BM_compatible: ClassVar[Literal[False]]
    def match(self, node): ...
    def transform(self, node, results): ...
