import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._tickfont import TickfontValidator
    from ._stream import StreamValidator
    from ._sortpaths import SortpathsValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._labelfont import LabelfontValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoveron import HoveronValidator
    from ._hoverinfo import HoverinfoValidator
    from ._domain import DomainValidator
    from ._dimensiondefaults import DimensiondefaultsValidator
    from ._dimensions import DimensionsValidator
    from ._countssrc import CountssrcValidator
    from ._counts import CountsValidator
    from ._bundlecolors import BundlecolorsValidator
    from ._arrangement import ArrangementValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._tickfont.TickfontValidator",
            "._stream.StreamValidator",
            "._sortpaths.SortpathsValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._labelfont.LabelfontValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoveron.HoveronValidator",
            "._hoverinfo.HoverinfoValidator",
            "._domain.DomainValidator",
            "._dimensiondefaults.DimensiondefaultsValidator",
            "._dimensions.DimensionsValidator",
            "._countssrc.CountssrcValidator",
            "._counts.CountsValidator",
            "._bundlecolors.BundlecolorsValidator",
            "._arrangement.ArrangementValidator",
        ],
    )
