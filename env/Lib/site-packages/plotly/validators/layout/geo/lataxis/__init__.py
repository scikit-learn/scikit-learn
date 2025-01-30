import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._tick0 import Tick0Validator
    from ._showgrid import ShowgridValidator
    from ._range import RangeValidator
    from ._gridwidth import GridwidthValidator
    from ._griddash import GriddashValidator
    from ._gridcolor import GridcolorValidator
    from ._dtick import DtickValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._tick0.Tick0Validator",
            "._showgrid.ShowgridValidator",
            "._range.RangeValidator",
            "._gridwidth.GridwidthValidator",
            "._griddash.GriddashValidator",
            "._gridcolor.GridcolorValidator",
            "._dtick.DtickValidator",
        ],
    )
