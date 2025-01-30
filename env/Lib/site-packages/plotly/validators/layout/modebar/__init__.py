import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._uirevision import UirevisionValidator
    from ._removesrc import RemovesrcValidator
    from ._remove import RemoveValidator
    from ._orientation import OrientationValidator
    from ._color import ColorValidator
    from ._bgcolor import BgcolorValidator
    from ._addsrc import AddsrcValidator
    from ._add import AddValidator
    from ._activecolor import ActivecolorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._uirevision.UirevisionValidator",
            "._removesrc.RemovesrcValidator",
            "._remove.RemoveValidator",
            "._orientation.OrientationValidator",
            "._color.ColorValidator",
            "._bgcolor.BgcolorValidator",
            "._addsrc.AddsrcValidator",
            "._add.AddValidator",
            "._activecolor.ActivecolorValidator",
        ],
    )
