import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zsrc import ZsrcValidator
    from ._z import ZValidator
    from ._ysrc import YsrcValidator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._x import XValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zsrc.ZsrcValidator",
            "._z.ZValidator",
            "._ysrc.YsrcValidator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._x.XValidator",
        ],
    )
