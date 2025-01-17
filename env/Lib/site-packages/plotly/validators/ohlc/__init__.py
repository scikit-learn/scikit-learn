import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zorder import ZorderValidator
    from ._yhoverformat import YhoverformatValidator
    from ._yaxis import YaxisValidator
    from ._xsrc import XsrcValidator
    from ._xperiodalignment import XperiodalignmentValidator
    from ._xperiod0 import Xperiod0Validator
    from ._xperiod import XperiodValidator
    from ._xhoverformat import XhoverformatValidator
    from ._xcalendar import XcalendarValidator
    from ._xaxis import XaxisValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._tickwidth import TickwidthValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._showlegend import ShowlegendValidator
    from ._selectedpoints import SelectedpointsValidator
    from ._opensrc import OpensrcValidator
    from ._open import OpenValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._lowsrc import LowsrcValidator
    from ._low import LowValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._increasing import IncreasingValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._highsrc import HighsrcValidator
    from ._high import HighValidator
    from ._decreasing import DecreasingValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._closesrc import ClosesrcValidator
    from ._close import CloseValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zorder.ZorderValidator",
            "._yhoverformat.YhoverformatValidator",
            "._yaxis.YaxisValidator",
            "._xsrc.XsrcValidator",
            "._xperiodalignment.XperiodalignmentValidator",
            "._xperiod0.Xperiod0Validator",
            "._xperiod.XperiodValidator",
            "._xhoverformat.XhoverformatValidator",
            "._xcalendar.XcalendarValidator",
            "._xaxis.XaxisValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._tickwidth.TickwidthValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._showlegend.ShowlegendValidator",
            "._selectedpoints.SelectedpointsValidator",
            "._opensrc.OpensrcValidator",
            "._open.OpenValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._lowsrc.LowsrcValidator",
            "._low.LowValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._increasing.IncreasingValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._highsrc.HighsrcValidator",
            "._high.HighValidator",
            "._decreasing.DecreasingValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._closesrc.ClosesrcValidator",
            "._close.CloseValidator",
        ],
    )
