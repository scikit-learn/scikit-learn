import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._weight import WeightValidator
    from ._variant import VariantValidator
    from ._textcase import TextcaseValidator
    from ._style import StyleValidator
    from ._size import SizeValidator
    from ._shadow import ShadowValidator
    from ._lineposition import LinepositionValidator
    from ._family import FamilyValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._weight.WeightValidator",
            "._variant.VariantValidator",
            "._textcase.TextcaseValidator",
            "._style.StyleValidator",
            "._size.SizeValidator",
            "._shadow.ShadowValidator",
            "._lineposition.LinepositionValidator",
            "._family.FamilyValidator",
            "._color.ColorValidator",
        ],
    )
