import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._widthsrc import WidthsrcValidator
    from ._width import WidthValidator
    from ._visible import VisibleValidator
    from ._unselected import UnselectedValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._thetaunit import ThetaunitValidator
    from ._thetasrc import ThetasrcValidator
    from ._theta0 import Theta0Validator
    from ._theta import ThetaValidator
    from ._textsrc import TextsrcValidator
    from ._text import TextValidator
    from ._subplot import SubplotValidator
    from ._stream import StreamValidator
    from ._showlegend import ShowlegendValidator
    from ._selectedpoints import SelectedpointsValidator
    from ._selected import SelectedValidator
    from ._rsrc import RsrcValidator
    from ._r0 import R0Validator
    from ._r import RValidator
    from ._opacity import OpacityValidator
    from ._offsetsrc import OffsetsrcValidator
    from ._offset import OffsetValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._marker import MarkerValidator
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
    from ._dtheta import DthetaValidator
    from ._dr import DrValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._basesrc import BasesrcValidator
    from ._base import BaseValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._widthsrc.WidthsrcValidator",
            "._width.WidthValidator",
            "._visible.VisibleValidator",
            "._unselected.UnselectedValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._thetaunit.ThetaunitValidator",
            "._thetasrc.ThetasrcValidator",
            "._theta0.Theta0Validator",
            "._theta.ThetaValidator",
            "._textsrc.TextsrcValidator",
            "._text.TextValidator",
            "._subplot.SubplotValidator",
            "._stream.StreamValidator",
            "._showlegend.ShowlegendValidator",
            "._selectedpoints.SelectedpointsValidator",
            "._selected.SelectedValidator",
            "._rsrc.RsrcValidator",
            "._r0.R0Validator",
            "._r.RValidator",
            "._opacity.OpacityValidator",
            "._offsetsrc.OffsetsrcValidator",
            "._offset.OffsetValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._marker.MarkerValidator",
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
            "._dtheta.DthetaValidator",
            "._dr.DrValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._basesrc.BasesrcValidator",
            "._base.BaseValidator",
        ],
    )
