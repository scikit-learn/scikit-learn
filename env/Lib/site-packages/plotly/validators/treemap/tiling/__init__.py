import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._squarifyratio import SquarifyratioValidator
    from ._pad import PadValidator
    from ._packing import PackingValidator
    from ._flip import FlipValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._squarifyratio.SquarifyratioValidator",
            "._pad.PadValidator",
            "._packing.PackingValidator",
            "._flip.FlipValidator",
        ],
    )
