import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._thickness import ThicknessValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._range import RangeValidator
    from ._name import NameValidator
    from ._line import LineValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._thickness.ThicknessValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._range.RangeValidator",
            "._name.NameValidator",
            "._line.LineValidator",
            "._color.ColorValidator",
        ],
    )
