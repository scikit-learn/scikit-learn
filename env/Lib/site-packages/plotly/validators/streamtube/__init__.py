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
    from ._wsrc import WsrcValidator
    from ._whoverformat import WhoverformatValidator
    from ._w import WValidator
    from ._vsrc import VsrcValidator
    from ._visible import VisibleValidator
    from ._vhoverformat import VhoverformatValidator
    from ._v import VValidator
    from ._usrc import UsrcValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._uhoverformat import UhoverformatValidator
    from ._u import UValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._starts import StartsValidator
    from ._sizeref import SizerefValidator
    from ._showscale import ShowscaleValidator
    from ._showlegend import ShowlegendValidator
    from ._scene import SceneValidator
    from ._reversescale import ReversescaleValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._maxdisplayed import MaxdisplayedValidator
    from ._lightposition import LightpositionValidator
    from ._lighting import LightingValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._colorscale import ColorscaleValidator
    from ._colorbar import ColorbarValidator
    from ._coloraxis import ColoraxisValidator
    from ._cmin import CminValidator
    from ._cmid import CmidValidator
    from ._cmax import CmaxValidator
    from ._cauto import CautoValidator
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
            "._wsrc.WsrcValidator",
            "._whoverformat.WhoverformatValidator",
            "._w.WValidator",
            "._vsrc.VsrcValidator",
            "._visible.VisibleValidator",
            "._vhoverformat.VhoverformatValidator",
            "._v.VValidator",
            "._usrc.UsrcValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._uhoverformat.UhoverformatValidator",
            "._u.UValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._starts.StartsValidator",
            "._sizeref.SizerefValidator",
            "._showscale.ShowscaleValidator",
            "._showlegend.ShowlegendValidator",
            "._scene.SceneValidator",
            "._reversescale.ReversescaleValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._maxdisplayed.MaxdisplayedValidator",
            "._lightposition.LightpositionValidator",
            "._lighting.LightingValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._colorscale.ColorscaleValidator",
            "._colorbar.ColorbarValidator",
            "._coloraxis.ColoraxisValidator",
            "._cmin.CminValidator",
            "._cmid.CmidValidator",
            "._cmax.CmaxValidator",
            "._cauto.CautoValidator",
            "._autocolorscale.AutocolorscaleValidator",
        ],
    )
