import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._domain import Domain
    from ._hoverlabel import Hoverlabel
    from ._insidetextfont import Insidetextfont
    from ._leaf import Leaf
    from ._legendgrouptitle import Legendgrouptitle
    from ._marker import Marker
    from ._outsidetextfont import Outsidetextfont
    from ._root import Root
    from ._stream import Stream
    from ._textfont import Textfont
    from . import hoverlabel
    from . import legendgrouptitle
    from . import marker
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".marker"],
        [
            "._domain.Domain",
            "._hoverlabel.Hoverlabel",
            "._insidetextfont.Insidetextfont",
            "._leaf.Leaf",
            "._legendgrouptitle.Legendgrouptitle",
            "._marker.Marker",
            "._outsidetextfont.Outsidetextfont",
            "._root.Root",
            "._stream.Stream",
            "._textfont.Textfont",
        ],
    )
