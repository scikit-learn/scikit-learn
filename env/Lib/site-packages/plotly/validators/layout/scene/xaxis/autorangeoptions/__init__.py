import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._minallowed import MinallowedValidator
    from ._maxallowed import MaxallowedValidator
    from ._includesrc import IncludesrcValidator
    from ._include import IncludeValidator
    from ._clipmin import ClipminValidator
    from ._clipmax import ClipmaxValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._minallowed.MinallowedValidator",
            "._maxallowed.MaxallowedValidator",
            "._includesrc.IncludesrcValidator",
            "._include.IncludeValidator",
            "._clipmin.ClipminValidator",
            "._clipmax.ClipmaxValidator",
        ],
    )
