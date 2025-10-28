import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._bounds import Bounds
    from ._center import Center
    from ._domain import Domain
    from ._layer import Layer
    from . import layer
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".layer"],
        ["._bounds.Bounds", "._center.Center", "._domain.Domain", "._layer.Layer"],
    )
