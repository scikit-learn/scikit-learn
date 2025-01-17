import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._valuesuffix import ValuesuffixValidator
    from ._valueformat import ValueformatValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._textfont import TextfontValidator
    from ._stream import StreamValidator
    from ._selectedpoints import SelectedpointsValidator
    from ._orientation import OrientationValidator
    from ._node import NodeValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._link import LinkValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legend import LegendValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfo import HoverinfoValidator
    from ._domain import DomainValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._arrangement import ArrangementValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._valuesuffix.ValuesuffixValidator",
            "._valueformat.ValueformatValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._textfont.TextfontValidator",
            "._stream.StreamValidator",
            "._selectedpoints.SelectedpointsValidator",
            "._orientation.OrientationValidator",
            "._node.NodeValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._link.LinkValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legend.LegendValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfo.HoverinfoValidator",
            "._domain.DomainValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._arrangement.ArrangementValidator",
        ],
    )
