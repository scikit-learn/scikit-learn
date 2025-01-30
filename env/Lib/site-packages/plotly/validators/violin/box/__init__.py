import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._width import WidthValidator
    from ._visible import VisibleValidator
    from ._line import LineValidator
    from ._fillcolor import FillcolorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._width.WidthValidator",
            "._visible.VisibleValidator",
            "._line.LineValidator",
            "._fillcolor.FillcolorValidator",
        ],
    )
