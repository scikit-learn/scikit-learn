import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._namelengthsrc import NamelengthsrcValidator
    from ._namelength import NamelengthValidator
    from ._font import FontValidator
    from ._bordercolorsrc import BordercolorsrcValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolorsrc import BgcolorsrcValidator
    from ._bgcolor import BgcolorValidator
    from ._alignsrc import AlignsrcValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._namelengthsrc.NamelengthsrcValidator",
            "._namelength.NamelengthValidator",
            "._font.FontValidator",
            "._bordercolorsrc.BordercolorsrcValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolorsrc.BgcolorsrcValidator",
            "._bgcolor.BgcolorValidator",
            "._alignsrc.AlignsrcValidator",
            "._align.AlignValidator",
        ],
    )
