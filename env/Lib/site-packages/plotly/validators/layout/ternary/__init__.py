import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._uirevision import UirevisionValidator
    from ._sum import SumValidator
    from ._domain import DomainValidator
    from ._caxis import CaxisValidator
    from ._bgcolor import BgcolorValidator
    from ._baxis import BaxisValidator
    from ._aaxis import AaxisValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._uirevision.UirevisionValidator",
            "._sum.SumValidator",
            "._domain.DomainValidator",
            "._caxis.CaxisValidator",
            "._bgcolor.BgcolorValidator",
            "._baxis.BaxisValidator",
            "._aaxis.AaxisValidator",
        ],
    )
