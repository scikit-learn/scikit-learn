import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._stepsrc import StepsrcValidator
    from ._step import StepValidator
    from ._sizesrc import SizesrcValidator
    from ._size import SizeValidator
    from ._opacitysrc import OpacitysrcValidator
    from ._opacity import OpacityValidator
    from ._maxzoom import MaxzoomValidator
    from ._enabled import EnabledValidator
    from ._colorsrc import ColorsrcValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._stepsrc.StepsrcValidator",
            "._step.StepValidator",
            "._sizesrc.SizesrcValidator",
            "._size.SizeValidator",
            "._opacitysrc.OpacitysrcValidator",
            "._opacity.OpacityValidator",
            "._maxzoom.MaxzoomValidator",
            "._enabled.EnabledValidator",
            "._colorsrc.ColorsrcValidator",
            "._color.ColorValidator",
        ],
    )
