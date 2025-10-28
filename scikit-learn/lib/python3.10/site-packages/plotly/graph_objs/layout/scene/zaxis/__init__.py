import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._autorangeoptions import Autorangeoptions
    from ._tickfont import Tickfont
    from ._tickformatstop import Tickformatstop
    from ._title import Title
    from . import title
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".title"],
        [
            "._autorangeoptions.Autorangeoptions",
            "._tickfont.Tickfont",
            "._tickformatstop.Tickformatstop",
            "._title.Title",
        ],
    )
