import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zoom import ZoomValidator
    from ._uirevision import UirevisionValidator
    from ._style import StyleValidator
    from ._pitch import PitchValidator
    from ._layerdefaults import LayerdefaultsValidator
    from ._layers import LayersValidator
    from ._domain import DomainValidator
    from ._center import CenterValidator
    from ._bounds import BoundsValidator
    from ._bearing import BearingValidator
    from ._accesstoken import AccesstokenValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zoom.ZoomValidator",
            "._uirevision.UirevisionValidator",
            "._style.StyleValidator",
            "._pitch.PitchValidator",
            "._layerdefaults.LayerdefaultsValidator",
            "._layers.LayersValidator",
            "._domain.DomainValidator",
            "._center.CenterValidator",
            "._bounds.BoundsValidator",
            "._bearing.BearingValidator",
            "._accesstoken.AccesstokenValidator",
        ],
    )
