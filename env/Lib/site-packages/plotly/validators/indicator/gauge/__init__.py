import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._threshold import ThresholdValidator
    from ._stepdefaults import StepdefaultsValidator
    from ._steps import StepsValidator
    from ._shape import ShapeValidator
    from ._borderwidth import BorderwidthValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._bar import BarValidator
    from ._axis import AxisValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._threshold.ThresholdValidator",
            "._stepdefaults.StepdefaultsValidator",
            "._steps.StepsValidator",
            "._shape.ShapeValidator",
            "._borderwidth.BorderwidthValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._bar.BarValidator",
            "._axis.AxisValidator",
        ],
    )
