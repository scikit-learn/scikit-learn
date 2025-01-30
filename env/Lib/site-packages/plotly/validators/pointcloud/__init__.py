import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._ysrc import YsrcValidator
    from ._yboundssrc import YboundssrcValidator
    from ._ybounds import YboundsValidator
    from ._yaxis import YaxisValidator
    from ._y import YValidator
    from ._xysrc import XysrcValidator
    from ._xy import XyValidator
    from ._xsrc import XsrcValidator
    from ._xboundssrc import XboundssrcValidator
    from ._xbounds import XboundsValidator
    from ._xaxis import XaxisValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._showlegend import ShowlegendValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._marker import MarkerValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._indicessrc import IndicessrcValidator
    from ._indices import IndicesValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._ysrc.YsrcValidator",
            "._yboundssrc.YboundssrcValidator",
            "._ybounds.YboundsValidator",
            "._yaxis.YaxisValidator",
            "._y.YValidator",
            "._xysrc.XysrcValidator",
            "._xy.XyValidator",
            "._xsrc.XsrcValidator",
            "._xboundssrc.XboundssrcValidator",
            "._xbounds.XboundsValidator",
            "._xaxis.XaxisValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._showlegend.ShowlegendValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._marker.MarkerValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._indicessrc.IndicessrcValidator",
            "._indices.IndicesValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
        ],
    )
