import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._valuessrc import ValuessrcValidator
    from ._values import ValuesValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._title import TitleValidator
    from ._texttemplatesrc import TexttemplatesrcValidator
    from ._texttemplate import TexttemplateValidator
    from ._textsrc import TextsrcValidator
    from ._textpositionsrc import TextpositionsrcValidator
    from ._textposition import TextpositionValidator
    from ._textinfo import TextinfoValidator
    from ._textfont import TextfontValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._sort import SortValidator
    from ._showlegend import ShowlegendValidator
    from ._scalegroup import ScalegroupValidator
    from ._rotation import RotationValidator
    from ._pullsrc import PullsrcValidator
    from ._pull import PullValidator
    from ._outsidetextfont import OutsidetextfontValidator
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
    from ._labelssrc import LabelssrcValidator
    from ._labels import LabelsValidator
    from ._label0 import Label0Validator
    from ._insidetextorientation import InsidetextorientationValidator
    from ._insidetextfont import InsidetextfontValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hovertextsrc import HovertextsrcValidator
    from ._hovertext import HovertextValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfosrc import HoverinfosrcValidator
    from ._hoverinfo import HoverinfoValidator
    from ._hole import HoleValidator
    from ._domain import DomainValidator
    from ._dlabel import DlabelValidator
    from ._direction import DirectionValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._automargin import AutomarginValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._valuessrc.ValuessrcValidator",
            "._values.ValuesValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._title.TitleValidator",
            "._texttemplatesrc.TexttemplatesrcValidator",
            "._texttemplate.TexttemplateValidator",
            "._textsrc.TextsrcValidator",
            "._textpositionsrc.TextpositionsrcValidator",
            "._textposition.TextpositionValidator",
            "._textinfo.TextinfoValidator",
            "._textfont.TextfontValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._sort.SortValidator",
            "._showlegend.ShowlegendValidator",
            "._scalegroup.ScalegroupValidator",
            "._rotation.RotationValidator",
            "._pullsrc.PullsrcValidator",
            "._pull.PullValidator",
            "._outsidetextfont.OutsidetextfontValidator",
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
            "._labelssrc.LabelssrcValidator",
            "._labels.LabelsValidator",
            "._label0.Label0Validator",
            "._insidetextorientation.InsidetextorientationValidator",
            "._insidetextfont.InsidetextfontValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hovertextsrc.HovertextsrcValidator",
            "._hovertext.HovertextValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfosrc.HoverinfosrcValidator",
            "._hoverinfo.HoverinfoValidator",
            "._hole.HoleValidator",
            "._domain.DomainValidator",
            "._dlabel.DlabelValidator",
            "._direction.DirectionValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._automargin.AutomarginValidator",
        ],
    )
