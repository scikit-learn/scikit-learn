import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._tickfont import Tickfont
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__, [], ["._tickfont.Tickfont"]
    )
