import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._colorbar import ColorBar
    from ._contours import Contours
    from ._legendgrouptitle import Legendgrouptitle
    from ._line import Line
    from ._stream import Stream
    from . import colorbar
    from . import contours
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".colorbar", ".contours", ".legendgrouptitle"],
        [
            "._colorbar.ColorBar",
            "._contours.Contours",
            "._legendgrouptitle.Legendgrouptitle",
            "._line.Line",
            "._stream.Stream",
        ],
    )
