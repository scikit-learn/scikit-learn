import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._autorangeoptions import Autorangeoptions
    from ._minor import Minor
    from ._rangebreak import Rangebreak
    from ._rangeselector import Rangeselector
    from ._rangeslider import Rangeslider
    from ._tickfont import Tickfont
    from ._tickformatstop import Tickformatstop
    from ._title import Title
    from ._unifiedhovertitle import Unifiedhovertitle
    from . import rangeselector
    from . import rangeslider
    from . import title
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".rangeselector", ".rangeslider", ".title"],
        [
            "._autorangeoptions.Autorangeoptions",
            "._minor.Minor",
            "._rangebreak.Rangebreak",
            "._rangeselector.Rangeselector",
            "._rangeslider.Rangeslider",
            "._tickfont.Tickfont",
            "._tickformatstop.Tickformatstop",
            "._title.Title",
            "._unifiedhovertitle.Unifiedhovertitle",
        ],
    )
