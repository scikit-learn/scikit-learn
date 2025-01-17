import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._decreasing import Decreasing
    from ._hoverlabel import Hoverlabel
    from ._increasing import Increasing
    from ._legendgrouptitle import Legendgrouptitle
    from ._line import Line
    from ._stream import Stream
    from . import decreasing
    from . import hoverlabel
    from . import increasing
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".decreasing", ".hoverlabel", ".increasing", ".legendgrouptitle"],
        [
            "._decreasing.Decreasing",
            "._hoverlabel.Hoverlabel",
            "._increasing.Increasing",
            "._legendgrouptitle.Legendgrouptitle",
            "._line.Line",
            "._stream.Stream",
        ],
    )
