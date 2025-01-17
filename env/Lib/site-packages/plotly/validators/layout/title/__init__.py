import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yref import YrefValidator
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xref import XrefValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._text import TextValidator
    from ._subtitle import SubtitleValidator
    from ._pad import PadValidator
    from ._font import FontValidator
    from ._automargin import AutomarginValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yref.YrefValidator",
            "._yanchor.YanchorValidator",
            "._y.YValidator",
            "._xref.XrefValidator",
            "._xanchor.XanchorValidator",
            "._x.XValidator",
            "._text.TextValidator",
            "._subtitle.SubtitleValidator",
            "._pad.PadValidator",
            "._font.FontValidator",
            "._automargin.AutomarginValidator",
        ],
    )
