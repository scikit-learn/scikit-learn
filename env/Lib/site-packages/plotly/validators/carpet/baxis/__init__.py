import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._type import TypeValidator
    from ._title import TitleValidator
    from ._tickvalssrc import TickvalssrcValidator
    from ._tickvals import TickvalsValidator
    from ._ticktextsrc import TicktextsrcValidator
    from ._ticktext import TicktextValidator
    from ._ticksuffix import TicksuffixValidator
    from ._tickprefix import TickprefixValidator
    from ._tickmode import TickmodeValidator
    from ._tickformatstopdefaults import TickformatstopdefaultsValidator
    from ._tickformatstops import TickformatstopsValidator
    from ._tickformat import TickformatValidator
    from ._tickfont import TickfontValidator
    from ._tickangle import TickangleValidator
    from ._tick0 import Tick0Validator
    from ._startlinewidth import StartlinewidthValidator
    from ._startlinecolor import StartlinecolorValidator
    from ._startline import StartlineValidator
    from ._smoothing import SmoothingValidator
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
    from ._minorgridwidth import MinorgridwidthValidator
    from ._minorgriddash import MinorgriddashValidator
    from ._minorgridcount import MinorgridcountValidator
    from ._minorgridcolor import MinorgridcolorValidator
    from ._minexponent import MinexponentValidator
    from ._linewidth import LinewidthValidator
    from ._linecolor import LinecolorValidator
    from ._labelsuffix import LabelsuffixValidator
    from ._labelprefix import LabelprefixValidator
    from ._labelpadding import LabelpaddingValidator
    from ._labelalias import LabelaliasValidator
    from ._gridwidth import GridwidthValidator
    from ._griddash import GriddashValidator
    from ._gridcolor import GridcolorValidator
    from ._fixedrange import FixedrangeValidator
    from ._exponentformat import ExponentformatValidator
    from ._endlinewidth import EndlinewidthValidator
    from ._endlinecolor import EndlinecolorValidator
    from ._endline import EndlineValidator
    from ._dtick import DtickValidator
    from ._color import ColorValidator
    from ._cheatertype import CheatertypeValidator
    from ._categoryorder import CategoryorderValidator
    from ._categoryarraysrc import CategoryarraysrcValidator
    from ._categoryarray import CategoryarrayValidator
    from ._autotypenumbers import AutotypenumbersValidator
    from ._autorange import AutorangeValidator
    from ._arraytick0 import Arraytick0Validator
    from ._arraydtick import ArraydtickValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._type.TypeValidator",
            "._title.TitleValidator",
            "._tickvalssrc.TickvalssrcValidator",
            "._tickvals.TickvalsValidator",
            "._ticktextsrc.TicktextsrcValidator",
            "._ticktext.TicktextValidator",
            "._ticksuffix.TicksuffixValidator",
            "._tickprefix.TickprefixValidator",
            "._tickmode.TickmodeValidator",
            "._tickformatstopdefaults.TickformatstopdefaultsValidator",
            "._tickformatstops.TickformatstopsValidator",
            "._tickformat.TickformatValidator",
            "._tickfont.TickfontValidator",
            "._tickangle.TickangleValidator",
            "._tick0.Tick0Validator",
            "._startlinewidth.StartlinewidthValidator",
            "._startlinecolor.StartlinecolorValidator",
            "._startline.StartlineValidator",
            "._smoothing.SmoothingValidator",
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
            "._minorgridwidth.MinorgridwidthValidator",
            "._minorgriddash.MinorgriddashValidator",
            "._minorgridcount.MinorgridcountValidator",
            "._minorgridcolor.MinorgridcolorValidator",
            "._minexponent.MinexponentValidator",
            "._linewidth.LinewidthValidator",
            "._linecolor.LinecolorValidator",
            "._labelsuffix.LabelsuffixValidator",
            "._labelprefix.LabelprefixValidator",
            "._labelpadding.LabelpaddingValidator",
            "._labelalias.LabelaliasValidator",
            "._gridwidth.GridwidthValidator",
            "._griddash.GriddashValidator",
            "._gridcolor.GridcolorValidator",
            "._fixedrange.FixedrangeValidator",
            "._exponentformat.ExponentformatValidator",
            "._endlinewidth.EndlinewidthValidator",
            "._endlinecolor.EndlinecolorValidator",
            "._endline.EndlineValidator",
            "._dtick.DtickValidator",
            "._color.ColorValidator",
            "._cheatertype.CheatertypeValidator",
            "._categoryorder.CategoryorderValidator",
            "._categoryarraysrc.CategoryarraysrcValidator",
            "._categoryarray.CategoryarrayValidator",
            "._autotypenumbers.AutotypenumbersValidator",
            "._autorange.AutorangeValidator",
            "._arraytick0.Arraytick0Validator",
            "._arraydtick.ArraydtickValidator",
        ],
    )
