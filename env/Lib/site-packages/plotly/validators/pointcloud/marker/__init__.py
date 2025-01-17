import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._sizemin import SizeminValidator
    from ._sizemax import SizemaxValidator
    from ._opacity import OpacityValidator
    from ._color import ColorValidator
    from ._border import BorderValidator
    from ._blend import BlendValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._sizemin.SizeminValidator",
            "._sizemax.SizemaxValidator",
            "._opacity.OpacityValidator",
            "._color.ColorValidator",
            "._border.BorderValidator",
            "._blend.BlendValidator",
        ],
    )
