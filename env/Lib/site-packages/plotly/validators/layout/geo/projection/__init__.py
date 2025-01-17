import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._type import TypeValidator
    from ._tilt import TiltValidator
    from ._scale import ScaleValidator
    from ._rotation import RotationValidator
    from ._parallels import ParallelsValidator
    from ._distance import DistanceValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._type.TypeValidator",
            "._tilt.TiltValidator",
            "._scale.ScaleValidator",
            "._rotation.RotationValidator",
            "._parallels.ParallelsValidator",
            "._distance.DistanceValidator",
        ],
    )
