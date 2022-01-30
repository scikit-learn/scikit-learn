# pyright: reportIncompleteStub = false
from typing import Any

# note: this is a lie to make type checkers happy (they special
# case property). cache_readonly uses attribute names similar to
# property (fget) but it does not provide fset and fdel.
cache_readonly = property

def __getattr__(name: str) -> Any: ...  # incomplete
