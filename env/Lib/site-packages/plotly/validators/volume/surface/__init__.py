import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._show import ShowValidator
    from ._pattern import PatternValidator
    from ._fill import FillValidator
    from ._count import CountValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._show.ShowValidator",
            "._pattern.PatternValidator",
            "._fill.FillValidator",
            "._count.CountValidator",
        ],
    )
