import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._currentvalue import Currentvalue
    from ._font import Font
    from ._pad import Pad
    from ._step import Step
    from ._transition import Transition
    from . import currentvalue
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".currentvalue"],
        [
            "._currentvalue.Currentvalue",
            "._font.Font",
            "._pad.Pad",
            "._step.Step",
            "._transition.Transition",
        ],
    )
