import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._domain import Domain
    from ._hoverlabel import Hoverlabel
    from ._insidetextfont import Insidetextfont
    from ._leaf import Leaf
    from ._legendgrouptitle import Legendgrouptitle
    from ._marker import Marker
    from ._outsidetextfont import Outsidetextfont
    from ._pathbar import Pathbar
    from ._root import Root
    from ._stream import Stream
    from ._textfont import Textfont
    from ._tiling import Tiling
    from . import hoverlabel
    from . import legendgrouptitle
    from . import marker
    from . import pathbar
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".marker", ".pathbar"],
        [
            "._domain.Domain",
            "._hoverlabel.Hoverlabel",
            "._insidetextfont.Insidetextfont",
            "._leaf.Leaf",
            "._legendgrouptitle.Legendgrouptitle",
            "._marker.Marker",
            "._outsidetextfont.Outsidetextfont",
            "._pathbar.Pathbar",
            "._root.Root",
            "._stream.Stream",
            "._textfont.Textfont",
            "._tiling.Tiling",
        ],
    )
