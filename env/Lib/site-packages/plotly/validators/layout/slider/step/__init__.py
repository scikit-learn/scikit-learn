import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._value import ValueValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._name import NameValidator
    from ._method import MethodValidator
    from ._label import LabelValidator
    from ._execute import ExecuteValidator
    from ._args import ArgsValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._value.ValueValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._name.NameValidator",
            "._method.MethodValidator",
            "._label.LabelValidator",
            "._execute.ExecuteValidator",
            "._args.ArgsValidator",
        ],
    )
