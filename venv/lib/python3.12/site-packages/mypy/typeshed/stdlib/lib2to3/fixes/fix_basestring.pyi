from typing import ClassVar, Literal

from .. import fixer_base

class FixBasestring(fixer_base.BaseFix):
    BM_compatible: ClassVar[Literal[True]]
    PATTERN: ClassVar[Literal["'basestring'"]]
    def transform(self, node, results): ...
