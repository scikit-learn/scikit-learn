import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._tickwidth import TickwidthValidator
    from ._tickvalssrc import TickvalssrcValidator
    from ._tickvals import TickvalsValidator
    from ._ticks import TicksValidator
    from ._tickmode import TickmodeValidator
    from ._ticklen import TicklenValidator
    from ._tickcolor import TickcolorValidator
    from ._tick0 import Tick0Validator
    from ._showgrid import ShowgridValidator
    from ._nticks import NticksValidator
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
            "._tickwidth.TickwidthValidator",
            "._tickvalssrc.TickvalssrcValidator",
            "._tickvals.TickvalsValidator",
            "._ticks.TicksValidator",
            "._tickmode.TickmodeValidator",
            "._ticklen.TicklenValidator",
            "._tickcolor.TickcolorValidator",
            "._tick0.Tick0Validator",
            "._showgrid.ShowgridValidator",
            "._nticks.NticksValidator",
            "._gridwidth.GridwidthValidator",
            "._griddash.GriddashValidator",
            "._gridcolor.GridcolorValidator",
            "._dtick.DtickValidator",
        ],
    )
