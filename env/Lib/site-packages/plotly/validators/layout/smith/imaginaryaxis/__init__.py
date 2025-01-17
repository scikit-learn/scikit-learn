import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._tickwidth import TickwidthValidator
    from ._tickvalssrc import TickvalssrcValidator
    from ._tickvals import TickvalsValidator
    from ._ticksuffix import TicksuffixValidator
    from ._ticks import TicksValidator
    from ._tickprefix import TickprefixValidator
    from ._ticklen import TicklenValidator
    from ._tickformat import TickformatValidator
    from ._tickfont import TickfontValidator
    from ._tickcolor import TickcolorValidator
    from ._showticksuffix import ShowticksuffixValidator
    from ._showtickprefix import ShowtickprefixValidator
    from ._showticklabels import ShowticklabelsValidator
    from ._showline import ShowlineValidator
    from ._showgrid import ShowgridValidator
    from ._linewidth import LinewidthValidator
    from ._linecolor import LinecolorValidator
    from ._layer import LayerValidator
    from ._labelalias import LabelaliasValidator
    from ._hoverformat import HoverformatValidator
    from ._gridwidth import GridwidthValidator
    from ._griddash import GriddashValidator
    from ._gridcolor import GridcolorValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._tickwidth.TickwidthValidator",
            "._tickvalssrc.TickvalssrcValidator",
            "._tickvals.TickvalsValidator",
            "._ticksuffix.TicksuffixValidator",
            "._ticks.TicksValidator",
            "._tickprefix.TickprefixValidator",
            "._ticklen.TicklenValidator",
            "._tickformat.TickformatValidator",
            "._tickfont.TickfontValidator",
            "._tickcolor.TickcolorValidator",
            "._showticksuffix.ShowticksuffixValidator",
            "._showtickprefix.ShowtickprefixValidator",
            "._showticklabels.ShowticklabelsValidator",
            "._showline.ShowlineValidator",
            "._showgrid.ShowgridValidator",
            "._linewidth.LinewidthValidator",
            "._linecolor.LinecolorValidator",
            "._layer.LayerValidator",
            "._labelalias.LabelaliasValidator",
            "._hoverformat.HoverformatValidator",
            "._gridwidth.GridwidthValidator",
            "._griddash.GriddashValidator",
            "._gridcolor.GridcolorValidator",
            "._color.ColorValidator",
        ],
    )
