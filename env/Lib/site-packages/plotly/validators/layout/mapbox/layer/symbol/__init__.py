import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._textposition import TextpositionValidator
    from ._textfont import TextfontValidator
    from ._text import TextValidator
    from ._placement import PlacementValidator
    from ._iconsize import IconsizeValidator
    from ._icon import IconValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._textposition.TextpositionValidator",
            "._textfont.TextfontValidator",
            "._text.TextValidator",
            "._placement.PlacementValidator",
            "._iconsize.IconsizeValidator",
            "._icon.IconValidator",
        ],
    )
