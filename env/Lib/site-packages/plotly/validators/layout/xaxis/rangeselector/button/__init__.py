import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._stepmode import StepmodeValidator
    from ._step import StepValidator
    from ._name import NameValidator
    from ._label import LabelValidator
    from ._count import CountValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._stepmode.StepmodeValidator",
            "._step.StepValidator",
            "._name.NameValidator",
            "._label.LabelValidator",
            "._count.CountValidator",
        ],
    )
