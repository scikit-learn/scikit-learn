import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._font import Font
    from ._pad import Pad
    from ._subtitle import Subtitle
    from . import subtitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__, [".subtitle"], ["._font.Font", "._pad.Pad", "._subtitle.Subtitle"]
    )
