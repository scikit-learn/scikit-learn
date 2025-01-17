import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._value import ValueValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._title import TitleValidator
    from ._stream import StreamValidator
    from ._number import NumberValidator
    from ._name import NameValidator
    from ._mode import ModeValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legend import LegendValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._gauge import GaugeValidator
    from ._domain import DomainValidator
    from ._delta import DeltaValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._value.ValueValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._title.TitleValidator",
            "._stream.StreamValidator",
            "._number.NumberValidator",
            "._name.NameValidator",
            "._mode.ModeValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legend.LegendValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._gauge.GaugeValidator",
            "._domain.DomainValidator",
            "._delta.DeltaValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._align.AlignValidator",
        ],
    )
