import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._connector import Connector
    from ._decreasing import Decreasing
    from ._hoverlabel import Hoverlabel
    from ._increasing import Increasing
    from ._insidetextfont import Insidetextfont
    from ._legendgrouptitle import Legendgrouptitle
    from ._outsidetextfont import Outsidetextfont
    from ._stream import Stream
    from ._textfont import Textfont
    from ._totals import Totals
    from . import connector
    from . import decreasing
    from . import hoverlabel
    from . import increasing
    from . import legendgrouptitle
    from . import totals
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [
            ".connector",
            ".decreasing",
            ".hoverlabel",
            ".increasing",
            ".legendgrouptitle",
            ".totals",
        ],
        [
            "._connector.Connector",
            "._decreasing.Decreasing",
            "._hoverlabel.Hoverlabel",
            "._increasing.Increasing",
            "._insidetextfont.Insidetextfont",
            "._legendgrouptitle.Legendgrouptitle",
            "._outsidetextfont.Outsidetextfont",
            "._stream.Stream",
            "._textfont.Textfont",
            "._totals.Totals",
        ],
    )
