import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._transition import TransitionValidator
    from ._tickwidth import TickwidthValidator
    from ._ticklen import TicklenValidator
    from ._tickcolor import TickcolorValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._stepdefaults import StepdefaultsValidator
    from ._steps import StepsValidator
    from ._pad import PadValidator
    from ._name import NameValidator
    from ._minorticklen import MinorticklenValidator
    from ._lenmode import LenmodeValidator
    from ._len import LenValidator
    from ._font import FontValidator
    from ._currentvalue import CurrentvalueValidator
    from ._borderwidth import BorderwidthValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._activebgcolor import ActivebgcolorValidator
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
            "._transition.TransitionValidator",
            "._tickwidth.TickwidthValidator",
            "._ticklen.TicklenValidator",
            "._tickcolor.TickcolorValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._stepdefaults.StepdefaultsValidator",
            "._steps.StepsValidator",
            "._pad.PadValidator",
            "._name.NameValidator",
            "._minorticklen.MinorticklenValidator",
            "._lenmode.LenmodeValidator",
            "._len.LenValidator",
            "._font.FontValidator",
            "._currentvalue.CurrentvalueValidator",
            "._borderwidth.BorderwidthValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._activebgcolor.ActivebgcolorValidator",
            "._active.ActiveValidator",
        ],
    )
