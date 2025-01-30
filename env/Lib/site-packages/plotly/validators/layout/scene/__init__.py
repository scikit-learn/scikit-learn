import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zaxis import ZaxisValidator
    from ._yaxis import YaxisValidator
    from ._xaxis import XaxisValidator
    from ._uirevision import UirevisionValidator
    from ._hovermode import HovermodeValidator
    from ._dragmode import DragmodeValidator
    from ._domain import DomainValidator
    from ._camera import CameraValidator
    from ._bgcolor import BgcolorValidator
    from ._aspectratio import AspectratioValidator
    from ._aspectmode import AspectmodeValidator
    from ._annotationdefaults import AnnotationdefaultsValidator
    from ._annotations import AnnotationsValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zaxis.ZaxisValidator",
            "._yaxis.YaxisValidator",
            "._xaxis.XaxisValidator",
            "._uirevision.UirevisionValidator",
            "._hovermode.HovermodeValidator",
            "._dragmode.DragmodeValidator",
            "._domain.DomainValidator",
            "._camera.CameraValidator",
            "._bgcolor.BgcolorValidator",
            "._aspectratio.AspectratioValidator",
            "._aspectmode.AspectmodeValidator",
            "._annotationdefaults.AnnotationdefaultsValidator",
            "._annotations.AnnotationsValidator",
        ],
    )
