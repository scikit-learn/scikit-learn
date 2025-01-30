import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._domain import Domain
    from ._hoverlabel import Hoverlabel
    from ._insidetextfont import Insidetextfont
    from ._legendgrouptitle import Legendgrouptitle
    from ._marker import Marker
    from ._outsidetextfont import Outsidetextfont
    from ._stream import Stream
    from ._textfont import Textfont
    from ._title import Title
    from . import hoverlabel
    from . import legendgrouptitle
    from . import marker
    from . import title
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".marker", ".title"],
        [
            "._domain.Domain",
            "._hoverlabel.Hoverlabel",
            "._insidetextfont.Insidetextfont",
            "._legendgrouptitle.Legendgrouptitle",
            "._marker.Marker",
            "._outsidetextfont.Outsidetextfont",
            "._stream.Stream",
            "._textfont.Textfont",
            "._title.Title",
        ],
    )
