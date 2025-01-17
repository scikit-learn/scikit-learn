import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._valueformat import ValueformatValidator
    from ._suffix import SuffixValidator
    from ._prefix import PrefixValidator
    from ._font import FontValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._valueformat.ValueformatValidator",
            "._suffix.SuffixValidator",
            "._prefix.PrefixValidator",
            "._font.FontValidator",
        ],
    )
