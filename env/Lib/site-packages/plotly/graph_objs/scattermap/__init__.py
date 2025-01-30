import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._cluster import Cluster
    from ._hoverlabel import Hoverlabel
    from ._legendgrouptitle import Legendgrouptitle
    from ._line import Line
    from ._marker import Marker
    from ._selected import Selected
    from ._stream import Stream
    from ._textfont import Textfont
    from ._unselected import Unselected
    from . import hoverlabel
    from . import legendgrouptitle
    from . import marker
    from . import selected
    from . import unselected
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".marker", ".selected", ".unselected"],
        [
            "._cluster.Cluster",
            "._hoverlabel.Hoverlabel",
            "._legendgrouptitle.Legendgrouptitle",
            "._line.Line",
            "._marker.Marker",
            "._selected.Selected",
            "._stream.Stream",
            "._textfont.Textfont",
            "._unselected.Unselected",
        ],
    )
