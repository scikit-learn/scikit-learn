import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._z import ZValidator
    from ._y import YValidator
    from ._x import XValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__, [], ["._z.ZValidator", "._y.YValidator", "._x.XValidator"]
    )
