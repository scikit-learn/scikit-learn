import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zsrc import ZsrcValidator
    from ._zhoverformat import ZhoverformatValidator
    from ._zcalendar import ZcalendarValidator
    from ._z import ZValidator
    from ._ysrc import YsrcValidator
    from ._yhoverformat import YhoverformatValidator
    from ._ycalendar import YcalendarValidator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._xhoverformat import XhoverformatValidator
    from ._xcalendar import XcalendarValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._texttemplatesrc import TexttemplatesrcValidator
    from ._texttemplate import TexttemplateValidator
    from ._textsrc import TextsrcValidator
    from ._textpositionsrc import TextpositionsrcValidator
    from ._textposition import TextpositionValidator
    from ._textfont import TextfontValidator
    from ._text import TextValidator
    from ._surfacecolor import SurfacecolorValidator
    from ._surfaceaxis import SurfaceaxisValidator
    from ._stream import StreamValidator
    from ._showlegend import ShowlegendValidator
    from ._scene import SceneValidator
    from ._projection import ProjectionValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._mode import ModeValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._marker import MarkerValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._error_z import Error_ZValidator
    from ._error_y import Error_YValidator
    from ._error_x import Error_XValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._connectgaps import ConnectgapsValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zsrc.ZsrcValidator",
            "._zhoverformat.ZhoverformatValidator",
            "._zcalendar.ZcalendarValidator",
            "._z.ZValidator",
            "._ysrc.YsrcValidator",
            "._yhoverformat.YhoverformatValidator",
            "._ycalendar.YcalendarValidator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._xhoverformat.XhoverformatValidator",
            "._xcalendar.XcalendarValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._texttemplatesrc.TexttemplatesrcValidator",
            "._texttemplate.TexttemplateValidator",
            "._textsrc.TextsrcValidator",
            "._textpositionsrc.TextpositionsrcValidator",
            "._textposition.TextpositionValidator",
            "._textfont.TextfontValidator",
            "._text.TextValidator",
            "._surfacecolor.SurfacecolorValidator",
            "._surfaceaxis.SurfaceaxisValidator",
            "._stream.StreamValidator",
            "._showlegend.ShowlegendValidator",
            "._scene.SceneValidator",
            "._projection.ProjectionValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._mode.ModeValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._marker.MarkerValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._error_z.Error_ZValidator",
            "._error_y.Error_YValidator",
            "._error_x.Error_XValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._connectgaps.ConnectgapsValidator",
        ],
    )
