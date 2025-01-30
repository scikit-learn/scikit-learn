import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._west import WestValidator
    from ._south import SouthValidator
    from ._north import NorthValidator
    from ._east import EastValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._west.WestValidator",
            "._south.SouthValidator",
            "._north.NorthValidator",
            "._east.EastValidator",
        ],
    )
