import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._center import Center
    from ._domain import Domain
    from ._lataxis import Lataxis
    from ._lonaxis import Lonaxis
    from ._projection import Projection
    from . import projection
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".projection"],
        [
            "._center.Center",
            "._domain.Domain",
            "._lataxis.Lataxis",
            "._lonaxis.Lonaxis",
            "._projection.Projection",
        ],
    )
