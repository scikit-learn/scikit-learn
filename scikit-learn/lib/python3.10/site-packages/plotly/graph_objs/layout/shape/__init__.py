import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._label import Label
    from ._legendgrouptitle import Legendgrouptitle
    from ._line import Line
    from . import label
    from . import legendgrouptitle
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".label", ".legendgrouptitle"],
        ["._label.Label", "._legendgrouptitle.Legendgrouptitle", "._line.Line"],
    )
