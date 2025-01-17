import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zsrc import ZsrcValidator
    from ._zhoverformat import ZhoverformatValidator
    from ._z import ZValidator
    from ._ysrc import YsrcValidator
    from ._yhoverformat import YhoverformatValidator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._xhoverformat import XhoverformatValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._valuesrc import ValuesrcValidator
    from ._valuehoverformat import ValuehoverformatValidator
    from ._value import ValueValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._surface import SurfaceValidator
    from ._stream import StreamValidator
    from ._spaceframe import SpaceframeValidator
    from ._slices import SlicesValidator
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
    from ._isomin import IsominValidator
    from ._isomax import IsomaxValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._flatshading import FlatshadingValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._contour import ContourValidator
    from ._colorscale import ColorscaleValidator
    from ._colorbar import ColorbarValidator
    from ._coloraxis import ColoraxisValidator
    from ._cmin import CminValidator
    from ._cmid import CmidValidator
    from ._cmax import CmaxValidator
    from ._cauto import CautoValidator
    from ._caps import CapsValidator
    from ._autocolorscale import AutocolorscaleValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zsrc.ZsrcValidator",
            "._zhoverformat.ZhoverformatValidator",
            "._z.ZValidator",
            "._ysrc.YsrcValidator",
            "._yhoverformat.YhoverformatValidator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._xhoverformat.XhoverformatValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._valuesrc.ValuesrcValidator",
            "._valuehoverformat.ValuehoverformatValidator",
            "._value.ValueValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._surface.SurfaceValidator",
            "._stream.StreamValidator",
            "._spaceframe.SpaceframeValidator",
            "._slices.SlicesValidator",
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
            "._isomin.IsominValidator",
            "._isomax.IsomaxValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._flatshading.FlatshadingValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._contour.ContourValidator",
            "._colorscale.ColorscaleValidator",
            "._colorbar.ColorbarValidator",
            "._coloraxis.ColoraxisValidator",
            "._cmin.CminValidator",
            "._cmid.CmidValidator",
            "._cmax.CmaxValidator",
            "._cauto.CautoValidator",
            "._caps.CapsValidator",
            "._autocolorscale.AutocolorscaleValidator",
        ],
    )
