import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._font import FontValidator
    from ._buttondefaults import ButtondefaultsValidator
    from ._buttons import ButtonsValidator
    from ._borderwidth import BorderwidthValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._activecolor import ActivecolorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yanchor.YanchorValidator",
            "._y.YValidator",
            "._xanchor.XanchorValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._font.FontValidator",
            "._buttondefaults.ButtondefaultsValidator",
            "._buttons.ButtonsValidator",
            "._borderwidth.BorderwidthValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._activecolor.ActivecolorValidator",
        ],
    )
