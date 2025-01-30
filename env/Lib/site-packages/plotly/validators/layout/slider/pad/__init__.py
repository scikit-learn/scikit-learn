import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._t import TValidator
    from ._r import RValidator
    from ._l import LValidator
    from ._b import BValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        ["._t.TValidator", "._r.RValidator", "._l.LValidator", "._b.BValidator"],
    )
