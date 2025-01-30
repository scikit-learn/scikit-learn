import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._xanchor import XanchorValidator
    from ._visible import VisibleValidator
    from ._suffix import SuffixValidator
    from ._prefix import PrefixValidator
    from ._offset import OffsetValidator
    from ._font import FontValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._xanchor.XanchorValidator",
            "._visible.VisibleValidator",
            "._suffix.SuffixValidator",
            "._prefix.PrefixValidator",
            "._offset.OffsetValidator",
            "._font.FontValidator",
        ],
    )
