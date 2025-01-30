import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._valuessrc import ValuessrcValidator
    from ._values import ValuesValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._name import NameValidator
    from ._label import LabelValidator
    from ._axis import AxisValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._valuessrc.ValuessrcValidator",
            "._values.ValuesValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._name.NameValidator",
            "._label.LabelValidator",
            "._axis.AxisValidator",
        ],
    )
