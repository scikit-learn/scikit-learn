import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._valueformat import ValueformatValidator
    from ._suffix import SuffixValidator
    from ._relative import RelativeValidator
    from ._reference import ReferenceValidator
    from ._prefix import PrefixValidator
    from ._position import PositionValidator
    from ._increasing import IncreasingValidator
    from ._font import FontValidator
    from ._decreasing import DecreasingValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._valueformat.ValueformatValidator",
            "._suffix.SuffixValidator",
            "._relative.RelativeValidator",
            "._reference.ReferenceValidator",
            "._prefix.PrefixValidator",
            "._position.PositionValidator",
            "._increasing.IncreasingValidator",
            "._font.FontValidator",
            "._decreasing.DecreasingValidator",
        ],
    )
