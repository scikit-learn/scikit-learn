from typing import ClassVar, Literal

from _typeshed import Incomplete

from .. import fixer_base

context: Incomplete

class FixInput(fixer_base.BaseFix):
    BM_compatible: ClassVar[Literal[True]]
    PATTERN: ClassVar[str]
    def transform(self, node, results): ...
