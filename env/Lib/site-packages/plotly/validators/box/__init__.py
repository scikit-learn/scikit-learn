import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zorder import ZorderValidator
    from ._ysrc import YsrcValidator
    from ._yperiodalignment import YperiodalignmentValidator
    from ._yperiod0 import Yperiod0Validator
    from ._yperiod import YperiodValidator
    from ._yhoverformat import YhoverformatValidator
    from ._ycalendar import YcalendarValidator
    from ._yaxis import YaxisValidator
    from ._y0 import Y0Validator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._xperiodalignment import XperiodalignmentValidator
    from ._xperiod0 import Xperiod0Validator
    from ._xperiod import XperiodValidator
    from ._xhoverformat import XhoverformatValidator
    from ._xcalendar import XcalendarValidator
    from ._xaxis import XaxisValidator
    from ._x0 import X0Validator
    from ._x import XValidator
    from ._width import WidthValidator
    from ._whiskerwidth import WhiskerwidthValidator
    from ._visible import VisibleValidator
    from ._upperfencesrc import UpperfencesrcValidator
    from ._upperfence import UpperfenceValidator
    from ._unselected import UnselectedValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._sizemode import SizemodeValidator
    from ._showwhiskers import ShowwhiskersValidator
    from ._showlegend import ShowlegendValidator
    from ._selectedpoints import SelectedpointsValidator
    from ._selected import SelectedValidator
    from ._sdsrc import SdsrcValidator
    from ._sdmultiple import SdmultipleValidator
    from ._sd import SdValidator
    from ._quartilemethod import QuartilemethodValidator
    from ._q3src import Q3SrcValidator
    from ._q3 import Q3Validator
    from ._q1src import Q1SrcValidator
    from ._q1 import Q1Validator
    from ._pointpos import PointposValidator
    from ._orientation import OrientationValidator
    from ._opacity import OpacityValidator
    from ._offsetgroup import OffsetgroupValidator
    from ._notchwidth import NotchwidthValidator
    from ._notchspansrc import NotchspansrcValidator
    from ._notchspan import NotchspanValidator
    from ._notched import NotchedValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._mediansrc import MediansrcValidator
    from ._median import MedianValidator
    from ._meansrc import MeansrcValidator
    from ._mean import MeanValidator
    from ._marker import MarkerValidator
    from ._lowerfencesrc import LowerfencesrcValidator
    from ._lowerfence import LowerfenceValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._jitter import JitterValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoveron import HoveronValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._fillcolor import FillcolorValidator
    from ._dy import DyValidator
    from ._dx import DxValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._boxpoints import BoxpointsValidator
    from ._boxmean import BoxmeanValidator
    from ._alignmentgroup import AlignmentgroupValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zorder.ZorderValidator",
            "._ysrc.YsrcValidator",
            "._yperiodalignment.YperiodalignmentValidator",
            "._yperiod0.Yperiod0Validator",
            "._yperiod.YperiodValidator",
            "._yhoverformat.YhoverformatValidator",
            "._ycalendar.YcalendarValidator",
            "._yaxis.YaxisValidator",
            "._y0.Y0Validator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._xperiodalignment.XperiodalignmentValidator",
            "._xperiod0.Xperiod0Validator",
            "._xperiod.XperiodValidator",
            "._xhoverformat.XhoverformatValidator",
            "._xcalendar.XcalendarValidator",
            "._xaxis.XaxisValidator",
            "._x0.X0Validator",
            "._x.XValidator",
            "._width.WidthValidator",
            "._whiskerwidth.WhiskerwidthValidator",
            "._visible.VisibleValidator",
            "._upperfencesrc.UpperfencesrcValidator",
            "._upperfence.UpperfenceValidator",
            "._unselected.UnselectedValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._sizemode.SizemodeValidator",
            "._showwhiskers.ShowwhiskersValidator",
            "._showlegend.ShowlegendValidator",
            "._selectedpoints.SelectedpointsValidator",
            "._selected.SelectedValidator",
            "._sdsrc.SdsrcValidator",
            "._sdmultiple.SdmultipleValidator",
            "._sd.SdValidator",
            "._quartilemethod.QuartilemethodValidator",
            "._q3src.Q3SrcValidator",
            "._q3.Q3Validator",
            "._q1src.Q1SrcValidator",
            "._q1.Q1Validator",
            "._pointpos.PointposValidator",
            "._orientation.OrientationValidator",
            "._opacity.OpacityValidator",
            "._offsetgroup.OffsetgroupValidator",
            "._notchwidth.NotchwidthValidator",
            "._notchspansrc.NotchspansrcValidator",
            "._notchspan.NotchspanValidator",
            "._notched.NotchedValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._mediansrc.MediansrcValidator",
            "._median.MedianValidator",
            "._meansrc.MeansrcValidator",
            "._mean.MeanValidator",
            "._marker.MarkerValidator",
            "._lowerfencesrc.LowerfencesrcValidator",
            "._lowerfence.LowerfenceValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._jitter.JitterValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoveron.HoveronValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._fillcolor.FillcolorValidator",
            "._dy.DyValidator",
            "._dx.DxValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._boxpoints.BoxpointsValidator",
            "._boxmean.BoxmeanValidator",
            "._alignmentgroup.AlignmentgroupValidator",
        ],
    )
