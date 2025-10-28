import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._font import Font
    from ._grouptitlefont import Grouptitlefont
    from ._title import Title
    from . import title
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".title"],
        ["._font.Font", "._grouptitlefont.Grouptitlefont", "._title.Title"],
    )
