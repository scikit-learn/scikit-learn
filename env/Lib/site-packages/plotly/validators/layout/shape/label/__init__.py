import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yanchor import YanchorValidator
    from ._xanchor import XanchorValidator
    from ._texttemplate import TexttemplateValidator
    from ._textposition import TextpositionValidator
    from ._textangle import TextangleValidator
    from ._text import TextValidator
    from ._padding import PaddingValidator
    from ._font import FontValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yanchor.YanchorValidator",
            "._xanchor.XanchorValidator",
            "._texttemplate.TexttemplateValidator",
            "._textposition.TextpositionValidator",
            "._textangle.TextangleValidator",
            "._text.TextValidator",
            "._padding.PaddingValidator",
            "._font.FontValidator",
        ],
    )
