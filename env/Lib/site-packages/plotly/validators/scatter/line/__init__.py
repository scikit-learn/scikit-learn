import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._width import WidthValidator
    from ._smoothing import SmoothingValidator
    from ._simplify import SimplifyValidator
    from ._shape import ShapeValidator
    from ._dash import DashValidator
    from ._color import ColorValidator
    from ._backoffsrc import BackoffsrcValidator
    from ._backoff import BackoffValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._width.WidthValidator",
            "._smoothing.SmoothingValidator",
            "._simplify.SimplifyValidator",
            "._shape.ShapeValidator",
            "._dash.DashValidator",
            "._color.ColorValidator",
            "._backoffsrc.BackoffsrcValidator",
            "._backoff.BackoffValidator",
        ],
    )
