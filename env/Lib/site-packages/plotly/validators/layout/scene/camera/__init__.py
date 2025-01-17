import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._up import UpValidator
    from ._projection import ProjectionValidator
    from ._eye import EyeValidator
    from ._center import CenterValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._up.UpValidator",
            "._projection.ProjectionValidator",
            "._eye.EyeValidator",
            "._center.CenterValidator",
        ],
    )
