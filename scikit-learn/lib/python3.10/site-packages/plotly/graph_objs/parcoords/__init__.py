import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._dimension import Dimension
    from ._domain import Domain
    from ._labelfont import Labelfont
    from ._legendgrouptitle import Legendgrouptitle
    from ._line import Line
    from ._rangefont import Rangefont
    from ._stream import Stream
    from ._tickfont import Tickfont
    from ._unselected import Unselected
    from . import legendgrouptitle
    from . import line
    from . import unselected
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".legendgrouptitle", ".line", ".unselected"],
        [
            "._dimension.Dimension",
            "._domain.Domain",
            "._labelfont.Labelfont",
            "._legendgrouptitle.Legendgrouptitle",
            "._line.Line",
            "._rangefont.Rangefont",
            "._stream.Stream",
            "._tickfont.Tickfont",
            "._unselected.Unselected",
        ],
    )
