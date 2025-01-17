import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._valuessrc import ValuessrcValidator
    from ._values import ValuesValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._tiling import TilingValidator
    from ._texttemplatesrc import TexttemplatesrcValidator
    from ._texttemplate import TexttemplateValidator
    from ._textsrc import TextsrcValidator
    from ._textposition import TextpositionValidator
    from ._textinfo import TextinfoValidator
    from ._textfont import TextfontValidator
    from ._text import TextValidator
    from ._stream import StreamValidator
    from ._sort import SortValidator
    from ._root import RootValidator
    from ._pathbar import PathbarValidator
    from ._parentssrc import ParentssrcValidator
    from ._parents import ParentsValidator
    from ._outsidetextfont import OutsidetextfontValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._maxdepth import MaxdepthValidator
    from ._marker import MarkerValidator
    from ._level import LevelValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legend import LegendValidator
    from ._leaf import LeafValidator
    from ._labelssrc import LabelssrcValidator
    from ._labels import LabelsValidator
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
    from ._domain import DomainValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._count import CountValidator
    from ._branchvalues import BranchvaluesValidator
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
            "._tiling.TilingValidator",
            "._texttemplatesrc.TexttemplatesrcValidator",
            "._texttemplate.TexttemplateValidator",
            "._textsrc.TextsrcValidator",
            "._textposition.TextpositionValidator",
            "._textinfo.TextinfoValidator",
            "._textfont.TextfontValidator",
            "._text.TextValidator",
            "._stream.StreamValidator",
            "._sort.SortValidator",
            "._root.RootValidator",
            "._pathbar.PathbarValidator",
            "._parentssrc.ParentssrcValidator",
            "._parents.ParentsValidator",
            "._outsidetextfont.OutsidetextfontValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._maxdepth.MaxdepthValidator",
            "._marker.MarkerValidator",
            "._level.LevelValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legend.LegendValidator",
            "._leaf.LeafValidator",
            "._labelssrc.LabelssrcValidator",
            "._labels.LabelsValidator",
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
            "._domain.DomainValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._count.CountValidator",
            "._branchvalues.BranchvaluesValidator",
        ],
    )
