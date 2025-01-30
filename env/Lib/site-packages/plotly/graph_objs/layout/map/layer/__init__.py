import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._circle import Circle
    from ._fill import Fill
    from ._line import Line
    from ._symbol import Symbol
    from . import symbol
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".symbol"],
        ["._circle.Circle", "._fill.Fill", "._line.Line", "._symbol.Symbol"],
    )
