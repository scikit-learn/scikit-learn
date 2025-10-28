import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._annotation import Annotation
    from ._aspectratio import Aspectratio
    from ._camera import Camera
    from ._domain import Domain
    from ._xaxis import XAxis
    from ._yaxis import YAxis
    from ._zaxis import ZAxis
    from . import annotation
    from . import camera
    from . import xaxis
    from . import yaxis
    from . import zaxis
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [".annotation", ".camera", ".xaxis", ".yaxis", ".zaxis"],
        [
            "._annotation.Annotation",
            "._aspectratio.Aspectratio",
            "._camera.Camera",
            "._domain.Domain",
            "._xaxis.XAxis",
            "._yaxis.YAxis",
            "._zaxis.ZAxis",
        ],
    )
