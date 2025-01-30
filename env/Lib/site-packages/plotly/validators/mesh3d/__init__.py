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
    from ._vertexcolorsrc import VertexcolorsrcValidator
    from ._vertexcolor import VertexcolorValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._showscale import ShowscaleValidator
    from ._showlegend import ShowlegendValidator
    from ._scene import SceneValidator
    from ._reversescale import ReversescaleValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._lightposition import LightpositionValidator
    from ._lighting import LightingValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._ksrc import KsrcValidator
    from ._k import KValidator
    from ._jsrc import JsrcValidator
    from ._j import JValidator
    from ._isrc import IsrcValidator
    from ._intensitysrc import IntensitysrcValidator
    from ._intensitymode import IntensitymodeValidator
    from ._intensity import IntensityValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._i import IValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._flatshading import FlatshadingValidator
    from ._facecolorsrc import FacecolorsrcValidator
    from ._facecolor import FacecolorValidator
    from ._delaunayaxis import DelaunayaxisValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._contour import ContourValidator
    from ._colorscale import ColorscaleValidator
    from ._colorbar import ColorbarValidator
    from ._coloraxis import ColoraxisValidator
    from ._color import ColorValidator
    from ._cmin import CminValidator
    from ._cmid import CmidValidator
    from ._cmax import CmaxValidator
    from ._cauto import CautoValidator
    from ._autocolorscale import AutocolorscaleValidator
    from ._alphahull import AlphahullValidator
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
            "._vertexcolorsrc.VertexcolorsrcValidator",
            "._vertexcolor.VertexcolorValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._showscale.ShowscaleValidator",
            "._showlegend.ShowlegendValidator",
            "._scene.SceneValidator",
            "._reversescale.ReversescaleValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._lightposition.LightpositionValidator",
            "._lighting.LightingValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._ksrc.KsrcValidator",
            "._k.KValidator",
            "._jsrc.JsrcValidator",
            "._j.JValidator",
            "._isrc.IsrcValidator",
            "._intensitysrc.IntensitysrcValidator",
            "._intensitymode.IntensitymodeValidator",
            "._intensity.IntensityValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._i.IValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._flatshading.FlatshadingValidator",
            "._facecolorsrc.FacecolorsrcValidator",
            "._facecolor.FacecolorValidator",
            "._delaunayaxis.DelaunayaxisValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._contour.ContourValidator",
            "._colorscale.ColorscaleValidator",
            "._colorbar.ColorbarValidator",
            "._coloraxis.ColoraxisValidator",
            "._color.ColorValidator",
            "._cmin.CminValidator",
            "._cmid.CmidValidator",
            "._cmax.CmaxValidator",
            "._cauto.CautoValidator",
            "._autocolorscale.AutocolorscaleValidator",
            "._alphahull.AlphahullValidator",
        ],
    )
