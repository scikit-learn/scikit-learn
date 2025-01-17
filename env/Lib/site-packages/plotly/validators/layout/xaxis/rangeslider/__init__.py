import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yaxis import YaxisValidator
    from ._visible import VisibleValidator
    from ._thickness import ThicknessValidator
    from ._range import RangeValidator
    from ._borderwidth import BorderwidthValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._autorange import AutorangeValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yaxis.YaxisValidator",
            "._visible.VisibleValidator",
            "._thickness.ThicknessValidator",
            "._range.RangeValidator",
            "._borderwidth.BorderwidthValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._autorange.AutorangeValidator",
        ],
    )
