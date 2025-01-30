import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._x import X
    from ._y import Y
    from ._z import Z
    from . import x
    from . import y
    from . import z
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__, [".x", ".y", ".z"], ["._x.X", "._y.Y", "._z.Z"]
    )
