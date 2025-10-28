import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._domain import Domain
    from ._hoverlabel import Hoverlabel
    from ._legendgrouptitle import Legendgrouptitle
    from ._link import Link
    from ._node import Node
    from ._stream import Stream
    from ._textfont import Textfont
    from . import hoverlabel
    from . import legendgrouptitle
    from . import link
    from . import node
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel", ".legendgrouptitle", ".link", ".node"],
        [
            "._domain.Domain",
            "._hoverlabel.Hoverlabel",
            "._legendgrouptitle.Legendgrouptitle",
            "._link.Link",
            "._node.Node",
            "._stream.Stream",
            "._textfont.Textfont",
        ],
    )
