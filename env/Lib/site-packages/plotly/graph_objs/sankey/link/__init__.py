import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._colorscale import Colorscale
    from ._hoverlabel import Hoverlabel
    from ._line import Line
    from . import hoverlabel
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".hoverlabel"],
        ["._colorscale.Colorscale", "._hoverlabel.Hoverlabel", "._line.Line"],
    )
