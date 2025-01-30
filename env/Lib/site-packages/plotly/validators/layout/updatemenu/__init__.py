import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._type import TypeValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._showactive import ShowactiveValidator
    from ._pad import PadValidator
    from ._name import NameValidator
    from ._font import FontValidator
    from ._direction import DirectionValidator
    from ._buttondefaults import ButtondefaultsValidator
    from ._buttons import ButtonsValidator
    from ._borderwidth import BorderwidthValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._active import ActiveValidator
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
            "._type.TypeValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._showactive.ShowactiveValidator",
            "._pad.PadValidator",
            "._name.NameValidator",
            "._font.FontValidator",
            "._direction.DirectionValidator",
            "._buttondefaults.ButtondefaultsValidator",
            "._buttons.ButtonsValidator",
            "._borderwidth.BorderwidthValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._active.ActiveValidator",
        ],
    )
