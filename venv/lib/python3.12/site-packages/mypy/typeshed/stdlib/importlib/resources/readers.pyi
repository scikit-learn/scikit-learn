# On py311+, things are actually defined here
# and re-exported from importlib.readers,
# but doing it this way leads to less code duplication for us

import sys
from collections.abc import Iterable, Iterator
from typing import TypeVar

if sys.version_info >= (3, 11):
    from importlib.readers import *

    _T = TypeVar("_T")

    def remove_duplicates(items: Iterable[_T]) -> Iterator[_T]: ...
