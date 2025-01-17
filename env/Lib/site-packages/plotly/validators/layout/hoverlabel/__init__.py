import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._namelength import NamelengthValidator
    from ._grouptitlefont import GrouptitlefontValidator
    from ._font import FontValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._namelength.NamelengthValidator",
            "._grouptitlefont.GrouptitlefontValidator",
            "._font.FontValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._align.AlignValidator",
        ],
    )
