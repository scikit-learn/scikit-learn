import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._specular import SpecularValidator
    from ._roughness import RoughnessValidator
    from ._fresnel import FresnelValidator
    from ._diffuse import DiffuseValidator
    from ._ambient import AmbientValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._specular.SpecularValidator",
            "._roughness.RoughnessValidator",
            "._fresnel.FresnelValidator",
            "._diffuse.DiffuseValidator",
            "._ambient.AmbientValidator",
        ],
    )
