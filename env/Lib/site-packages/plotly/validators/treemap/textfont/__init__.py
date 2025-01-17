import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._weightsrc import WeightsrcValidator
    from ._weight import WeightValidator
    from ._variantsrc import VariantsrcValidator
    from ._variant import VariantValidator
    from ._textcasesrc import TextcasesrcValidator
    from ._textcase import TextcaseValidator
    from ._stylesrc import StylesrcValidator
    from ._style import StyleValidator
    from ._sizesrc import SizesrcValidator
    from ._size import SizeValidator
    from ._shadowsrc import ShadowsrcValidator
    from ._shadow import ShadowValidator
    from ._linepositionsrc import LinepositionsrcValidator
    from ._lineposition import LinepositionValidator
    from ._familysrc import FamilysrcValidator
    from ._family import FamilyValidator
    from ._colorsrc import ColorsrcValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._weightsrc.WeightsrcValidator",
            "._weight.WeightValidator",
            "._variantsrc.VariantsrcValidator",
            "._variant.VariantValidator",
            "._textcasesrc.TextcasesrcValidator",
            "._textcase.TextcaseValidator",
            "._stylesrc.StylesrcValidator",
            "._style.StyleValidator",
            "._sizesrc.SizesrcValidator",
            "._size.SizeValidator",
            "._shadowsrc.ShadowsrcValidator",
            "._shadow.ShadowValidator",
            "._linepositionsrc.LinepositionsrcValidator",
            "._lineposition.LinepositionValidator",
            "._familysrc.FamilysrcValidator",
            "._family.FamilyValidator",
            "._colorsrc.ColorsrcValidator",
            "._color.ColorValidator",
        ],
    )
