import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._button import Button
    from ._font import Font
    from ._pad import Pad
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__, [], ["._button.Button", "._font.Font", "._pad.Pad"]
    )
