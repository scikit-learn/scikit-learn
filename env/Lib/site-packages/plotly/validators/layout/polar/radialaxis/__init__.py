import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._type import TypeValidator
    from ._title import TitleValidator
    from ._tickwidth import TickwidthValidator
    from ._tickvalssrc import TickvalssrcValidator
    from ._tickvals import TickvalsValidator
    from ._ticktextsrc import TicktextsrcValidator
    from ._ticktext import TicktextValidator
    from ._ticksuffix import TicksuffixValidator
    from ._ticks import TicksValidator
    from ._tickprefix import TickprefixValidator
    from ._tickmode import TickmodeValidator
    from ._ticklen import TicklenValidator
    from ._ticklabelstep import TicklabelstepValidator
    from ._tickformatstopdefaults import TickformatstopdefaultsValidator
    from ._tickformatstops import TickformatstopsValidator
    from ._tickformat import TickformatValidator
    from ._tickfont import TickfontValidator
    from ._tickcolor import TickcolorValidator
    from ._tickangle import TickangleValidator
    from ._tick0 import Tick0Validator
    from ._side import SideValidator
    from ._showticksuffix import ShowticksuffixValidator
    from ._showtickprefix import ShowtickprefixValidator
    from ._showticklabels import ShowticklabelsValidator
    from ._showline import ShowlineValidator
    from ._showgrid import ShowgridValidator
    from ._showexponent import ShowexponentValidator
    from ._separatethousands import SeparatethousandsValidator
    from ._rangemode import RangemodeValidator
    from ._range import RangeValidator
    from ._nticks import NticksValidator
    from ._minexponent import MinexponentValidator
    from ._minallowed import MinallowedValidator
    from ._maxallowed import MaxallowedValidator
    from ._linewidth import LinewidthValidator
    from ._linecolor import LinecolorValidator
    from ._layer import LayerValidator
    from ._labelalias import LabelaliasValidator
    from ._hoverformat import HoverformatValidator
    from ._gridwidth import GridwidthValidator
    from ._griddash import GriddashValidator
    from ._gridcolor import GridcolorValidator
    from ._exponentformat import ExponentformatValidator
    from ._dtick import DtickValidator
    from ._color import ColorValidator
    from ._categoryorder import CategoryorderValidator
    from ._categoryarraysrc import CategoryarraysrcValidator
    from ._categoryarray import CategoryarrayValidator
    from ._calendar import CalendarValidator
    from ._autotypenumbers import AutotypenumbersValidator
    from ._autotickangles import AutotickanglesValidator
    from ._autorangeoptions import AutorangeoptionsValidator
    from ._autorange import AutorangeValidator
    from ._angle import AngleValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._type.TypeValidator",
            "._title.TitleValidator",
            "._tickwidth.TickwidthValidator",
            "._tickvalssrc.TickvalssrcValidator",
            "._tickvals.TickvalsValidator",
            "._ticktextsrc.TicktextsrcValidator",
            "._ticktext.TicktextValidator",
            "._ticksuffix.TicksuffixValidator",
            "._ticks.TicksValidator",
            "._tickprefix.TickprefixValidator",
            "._tickmode.TickmodeValidator",
            "._ticklen.TicklenValidator",
            "._ticklabelstep.TicklabelstepValidator",
            "._tickformatstopdefaults.TickformatstopdefaultsValidator",
            "._tickformatstops.TickformatstopsValidator",
            "._tickformat.TickformatValidator",
            "._tickfont.TickfontValidator",
            "._tickcolor.TickcolorValidator",
            "._tickangle.TickangleValidator",
            "._tick0.Tick0Validator",
            "._side.SideValidator",
            "._showticksuffix.ShowticksuffixValidator",
            "._showtickprefix.ShowtickprefixValidator",
            "._showticklabels.ShowticklabelsValidator",
            "._showline.ShowlineValidator",
            "._showgrid.ShowgridValidator",
            "._showexponent.ShowexponentValidator",
            "._separatethousands.SeparatethousandsValidator",
            "._rangemode.RangemodeValidator",
            "._range.RangeValidator",
            "._nticks.NticksValidator",
            "._minexponent.MinexponentValidator",
            "._minallowed.MinallowedValidator",
            "._maxallowed.MaxallowedValidator",
            "._linewidth.LinewidthValidator",
            "._linecolor.LinecolorValidator",
            "._layer.LayerValidator",
            "._labelalias.LabelaliasValidator",
            "._hoverformat.HoverformatValidator",
            "._gridwidth.GridwidthValidator",
            "._griddash.GriddashValidator",
            "._gridcolor.GridcolorValidator",
            "._exponentformat.ExponentformatValidator",
            "._dtick.DtickValidator",
            "._color.ColorValidator",
            "._categoryorder.CategoryorderValidator",
            "._categoryarraysrc.CategoryarraysrcValidator",
            "._categoryarray.CategoryarrayValidator",
            "._calendar.CalendarValidator",
            "._autotypenumbers.AutotypenumbersValidator",
            "._autotickangles.AutotickanglesValidator",
            "._autorangeoptions.AutorangeoptionsValidator",
            "._autorange.AutorangeValidator",
            "._angle.AngleValidator",
        ],
    )
