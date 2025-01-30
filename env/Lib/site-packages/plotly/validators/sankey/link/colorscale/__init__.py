import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._templateitemname import TemplateitemnameValidator
    from ._name import NameValidator
    from ._label import LabelValidator
    from ._colorscale import ColorscaleValidator
    from ._cmin import CminValidator
    from ._cmax import CmaxValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._templateitemname.TemplateitemnameValidator",
            "._name.NameValidator",
            "._label.LabelValidator",
            "._colorscale.ColorscaleValidator",
            "._cmin.CminValidator",
            "._cmax.CmaxValidator",
        ],
    )
