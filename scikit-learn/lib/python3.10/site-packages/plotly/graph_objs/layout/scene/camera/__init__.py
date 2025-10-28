import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._center import Center
    from ._eye import Eye
    from ._projection import Projection
    from ._up import Up
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        ["._center.Center", "._eye.Eye", "._projection.Projection", "._up.Up"],
    )
