import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._colorbar import ColorBar
    from ._contours import Contours
    from ._hoverlabel import Hoverlabel
    from ._legendgrouptitle import Legendgrouptitle
    from ._lighting import Lighting
    from ._lightposition import Lightposition
    from ._stream import Stream
    from . import colorbar
    from . import contours
    from . import hoverlabel
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".colorbar", ".contours", ".hoverlabel", ".legendgrouptitle"],
        [
            "._colorbar.ColorBar",
            "._contours.Contours",
            "._hoverlabel.Hoverlabel",
            "._legendgrouptitle.Legendgrouptitle",
            "._lighting.Lighting",
            "._lightposition.Lightposition",
            "._stream.Stream",
        ],
    )
