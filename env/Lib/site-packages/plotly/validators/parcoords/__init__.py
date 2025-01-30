import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._unselected import UnselectedValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._tickfont import TickfontValidator
    from ._stream import StreamValidator
    from ._rangefont import RangefontValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legend import LegendValidator
    from ._labelside import LabelsideValidator
    from ._labelfont import LabelfontValidator
    from ._labelangle import LabelangleValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._domain import DomainValidator
    from ._dimensiondefaults import DimensiondefaultsValidator
    from ._dimensions import DimensionsValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._unselected.UnselectedValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._tickfont.TickfontValidator",
            "._stream.StreamValidator",
            "._rangefont.RangefontValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legend.LegendValidator",
            "._labelside.LabelsideValidator",
            "._labelfont.LabelfontValidator",
            "._labelangle.LabelangleValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._domain.DomainValidator",
            "._dimensiondefaults.DimensiondefaultsValidator",
            "._dimensions.DimensionsValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
        ],
    )
