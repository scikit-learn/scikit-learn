import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._symbolsrc import SymbolsrcValidator
    from ._symbol import SymbolValidator
    from ._sizesrc import SizesrcValidator
    from ._sizeref import SizerefValidator
    from ._sizemode import SizemodeValidator
    from ._sizemin import SizeminValidator
    from ._size import SizeValidator
    from ._showscale import ShowscaleValidator
    from ._reversescale import ReversescaleValidator
    from ._opacitysrc import OpacitysrcValidator
    from ._opacity import OpacityValidator
    from ._line import LineValidator
    from ._colorsrc import ColorsrcValidator
    from ._colorscale import ColorscaleValidator
    from ._colorbar import ColorbarValidator
    from ._coloraxis import ColoraxisValidator
    from ._color import ColorValidator
    from ._cmin import CminValidator
    from ._cmid import CmidValidator
    from ._cmax import CmaxValidator
    from ._cauto import CautoValidator
    from ._autocolorscale import AutocolorscaleValidator
    from ._anglesrc import AnglesrcValidator
    from ._angle import AngleValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._symbolsrc.SymbolsrcValidator",
            "._symbol.SymbolValidator",
            "._sizesrc.SizesrcValidator",
            "._sizeref.SizerefValidator",
            "._sizemode.SizemodeValidator",
            "._sizemin.SizeminValidator",
            "._size.SizeValidator",
            "._showscale.ShowscaleValidator",
            "._reversescale.ReversescaleValidator",
            "._opacitysrc.OpacitysrcValidator",
            "._opacity.OpacityValidator",
            "._line.LineValidator",
            "._colorsrc.ColorsrcValidator",
            "._colorscale.ColorscaleValidator",
            "._colorbar.ColorbarValidator",
            "._coloraxis.ColoraxisValidator",
            "._color.ColorValidator",
            "._cmin.CminValidator",
            "._cmid.CmidValidator",
            "._cmax.CmaxValidator",
            "._cauto.CautoValidator",
            "._autocolorscale.AutocolorscaleValidator",
            "._anglesrc.AnglesrcValidator",
            "._angle.AngleValidator",
        ],
    )
