import sys
from typing import List, Tuple

if sys.platform == "win32":

    ActionText: List[Tuple[str, str, str | None]]
    UIText: List[Tuple[str, str | None]]

    tables: List[str]
