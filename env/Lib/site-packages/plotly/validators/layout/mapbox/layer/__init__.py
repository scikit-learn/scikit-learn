import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._type import TypeValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._symbol import SymbolValidator
    from ._sourcetype import SourcetypeValidator
    from ._sourcelayer import SourcelayerValidator
    from ._sourceattribution import SourceattributionValidator
    from ._source import SourceValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._minzoom import MinzoomValidator
    from ._maxzoom import MaxzoomValidator
    from ._line import LineValidator
    from ._fill import FillValidator
    from ._coordinates import CoordinatesValidator
    from ._color import ColorValidator
    from ._circle import CircleValidator
    from ._below import BelowValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._type.TypeValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._symbol.SymbolValidator",
            "._sourcetype.SourcetypeValidator",
            "._sourcelayer.SourcelayerValidator",
            "._sourceattribution.SourceattributionValidator",
            "._source.SourceValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._minzoom.MinzoomValidator",
            "._maxzoom.MaxzoomValidator",
            "._line.LineValidator",
            "._fill.FillValidator",
            "._coordinates.CoordinatesValidator",
            "._color.ColorValidator",
            "._circle.CircleValidator",
            "._below.BelowValidator",
        ],
    )
