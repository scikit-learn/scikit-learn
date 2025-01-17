import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._uirevision import UirevisionValidator
    from ._sector import SectorValidator
    from ._radialaxis import RadialaxisValidator
    from ._hole import HoleValidator
    from ._gridshape import GridshapeValidator
    from ._domain import DomainValidator
    from ._bgcolor import BgcolorValidator
    from ._barmode import BarmodeValidator
    from ._bargap import BargapValidator
    from ._angularaxis import AngularaxisValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._uirevision.UirevisionValidator",
            "._sector.SectorValidator",
            "._radialaxis.RadialaxisValidator",
            "._hole.HoleValidator",
            "._gridshape.GridshapeValidator",
            "._domain.DomainValidator",
            "._bgcolor.BgcolorValidator",
            "._barmode.BarmodeValidator",
            "._bargap.BargapValidator",
            "._angularaxis.AngularaxisValidator",
        ],
    )
