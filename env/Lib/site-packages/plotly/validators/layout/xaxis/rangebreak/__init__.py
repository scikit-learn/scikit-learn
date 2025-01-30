import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._values import ValuesValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._pattern import PatternValidator
    from ._name import NameValidator
    from ._enabled import EnabledValidator
    from ._dvalue import DvalueValidator
    from ._bounds import BoundsValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._values.ValuesValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._pattern.PatternValidator",
            "._name.NameValidator",
            "._enabled.EnabledValidator",
            "._dvalue.DvalueValidator",
            "._bounds.BoundsValidator",
        ],
    )
