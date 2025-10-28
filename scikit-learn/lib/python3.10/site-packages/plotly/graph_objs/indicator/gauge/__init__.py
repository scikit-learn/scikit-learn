import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._axis import Axis
    from ._bar import Bar
    from ._step import Step
    from ._threshold import Threshold
    from . import axis
    from . import bar
    from . import step
    from . import threshold
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".axis", ".bar", ".step", ".threshold"],
        ["._axis.Axis", "._bar.Bar", "._step.Step", "._threshold.Threshold"],
    )
