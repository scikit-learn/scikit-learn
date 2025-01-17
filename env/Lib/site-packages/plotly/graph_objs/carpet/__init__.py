import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._aaxis import Aaxis
    from ._baxis import Baxis
    from ._font import Font
    from ._legendgrouptitle import Legendgrouptitle
    from ._stream import Stream
    from . import aaxis
    from . import baxis
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".aaxis", ".baxis", ".legendgrouptitle"],
        [
            "._aaxis.Aaxis",
            "._baxis.Baxis",
            "._font.Font",
            "._legendgrouptitle.Legendgrouptitle",
            "._stream.Stream",
        ],
    )
