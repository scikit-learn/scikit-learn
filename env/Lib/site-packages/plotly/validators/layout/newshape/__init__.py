import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._showlegend import ShowlegendValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._layer import LayerValidator
    from ._label import LabelValidator
    from ._fillrule import FillruleValidator
    from ._fillcolor import FillcolorValidator
    from ._drawdirection import DrawdirectionValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._showlegend.ShowlegendValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._layer.LayerValidator",
            "._label.LabelValidator",
            "._fillrule.FillruleValidator",
            "._fillcolor.FillcolorValidator",
            "._drawdirection.DrawdirectionValidator",
        ],
    )
