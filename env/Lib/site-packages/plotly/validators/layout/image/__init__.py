import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yref import YrefValidator
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xref import XrefValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._source import SourceValidator
    from ._sizing import SizingValidator
    from ._sizey import SizeyValidator
    from ._sizex import SizexValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._layer import LayerValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yref.YrefValidator",
            "._yanchor.YanchorValidator",
            "._y.YValidator",
            "._xref.XrefValidator",
            "._xanchor.XanchorValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._source.SourceValidator",
            "._sizing.SizingValidator",
            "._sizey.SizeyValidator",
            "._sizex.SizexValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._layer.LayerValidator",
        ],
    )
