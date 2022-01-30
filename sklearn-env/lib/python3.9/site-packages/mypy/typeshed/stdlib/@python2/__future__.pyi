import sys
from typing import List

class _Feature:
    def __init__(self, optionalRelease: sys._version_info, mandatoryRelease: sys._version_info, compiler_flag: int) -> None: ...
    def getOptionalRelease(self) -> sys._version_info: ...
    def getMandatoryRelease(self) -> sys._version_info: ...
    compiler_flag: int

absolute_import: _Feature
division: _Feature
generators: _Feature
nested_scopes: _Feature
print_function: _Feature
unicode_literals: _Feature
with_statement: _Feature
all_feature_names: List[str]  # undocumented
