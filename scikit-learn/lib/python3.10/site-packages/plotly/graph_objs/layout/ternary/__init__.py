import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._aaxis import Aaxis
    from ._baxis import Baxis
    from ._caxis import Caxis
    from ._domain import Domain
    from . import aaxis
    from . import baxis
    from . import caxis
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".aaxis", ".baxis", ".caxis"],
        ["._aaxis.Aaxis", "._baxis.Baxis", "._caxis.Caxis", "._domain.Domain"],
    )
