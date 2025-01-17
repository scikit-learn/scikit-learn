import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._colorbar import ColorBar
    from ._hoverlabel import Hoverlabel
    from ._legendgrouptitle import Legendgrouptitle
    from ._lighting import Lighting
    from ._lightposition import Lightposition
    from ._stream import Stream
    from . import colorbar
    from . import hoverlabel
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".colorbar", ".hoverlabel", ".legendgrouptitle"],
        [
            "._colorbar.ColorBar",
            "._hoverlabel.Hoverlabel",
            "._legendgrouptitle.Legendgrouptitle",
            "._lighting.Lighting",
            "._lightposition.Lightposition",
            "._stream.Stream",
        ],
    )
