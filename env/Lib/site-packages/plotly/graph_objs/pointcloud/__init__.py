import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._hoverlabel import Hoverlabel
    from ._legendgrouptitle import Legendgrouptitle
    from ._marker import Marker
    from ._stream import Stream
    from . import hoverlabel
    from . import legendgrouptitle
    from . import marker
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".marker"],
        [
            "._hoverlabel.Hoverlabel",
            "._legendgrouptitle.Legendgrouptitle",
            "._marker.Marker",
            "._stream.Stream",
        ],
    )
