import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._showscale import ShowscaleValidator
    from ._reversescale import ReversescaleValidator
    from ._pattern import PatternValidator
    from ._pad import PadValidator
    from ._line import LineValidator
    from ._depthfade import DepthfadeValidator
    from ._cornerradius import CornerradiusValidator
    from ._colorssrc import ColorssrcValidator
    from ._colorscale import ColorscaleValidator
    from ._colors import ColorsValidator
    from ._colorbar import ColorbarValidator
    from ._coloraxis import ColoraxisValidator
    from ._cmin import CminValidator
    from ._cmid import CmidValidator
    from ._cmax import CmaxValidator
    from ._cauto import CautoValidator
    from ._autocolorscale import AutocolorscaleValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._showscale.ShowscaleValidator",
            "._reversescale.ReversescaleValidator",
            "._pattern.PatternValidator",
            "._pad.PadValidator",
            "._line.LineValidator",
            "._depthfade.DepthfadeValidator",
            "._cornerradius.CornerradiusValidator",
            "._colorssrc.ColorssrcValidator",
            "._colorscale.ColorscaleValidator",
            "._colors.ColorsValidator",
            "._colorbar.ColorbarValidator",
            "._coloraxis.ColoraxisValidator",
            "._cmin.CminValidator",
            "._cmid.CmidValidator",
            "._cmax.CmaxValidator",
            "._cauto.CautoValidator",
            "._autocolorscale.AutocolorscaleValidator",
        ],
    )
