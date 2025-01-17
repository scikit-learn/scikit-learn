import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yref import YrefValidator
    from ._y1 import Y1Validator
    from ._y0 import Y0Validator
    from ._xref import XrefValidator
    from ._x1 import X1Validator
    from ._x0 import X0Validator
    from ._type import TypeValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._path import PathValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._line import LineValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yref.YrefValidator",
            "._y1.Y1Validator",
            "._y0.Y0Validator",
            "._xref.XrefValidator",
            "._x1.X1Validator",
            "._x0.X0Validator",
            "._type.TypeValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._path.PathValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._line.LineValidator",
        ],
    )
