import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._domain import Domain
    from ._imaginaryaxis import Imaginaryaxis
    from ._realaxis import Realaxis
    from . import imaginaryaxis
    from . import realaxis
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".imaginaryaxis", ".realaxis"],
        ["._domain.Domain", "._imaginaryaxis.Imaginaryaxis", "._realaxis.Realaxis"],
    )
