import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._valuessrc import ValuessrcValidator
    from ._values import ValuesValidator
    from ._suffixsrc import SuffixsrcValidator
    from ._suffix import SuffixValidator
    from ._prefixsrc import PrefixsrcValidator
    from ._prefix import PrefixValidator
    from ._line import LineValidator
    from ._height import HeightValidator
    from ._formatsrc import FormatsrcValidator
    from ._format import FormatValidator
    from ._font import FontValidator
    from ._fill import FillValidator
    from ._alignsrc import AlignsrcValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._valuessrc.ValuessrcValidator",
            "._values.ValuesValidator",
            "._suffixsrc.SuffixsrcValidator",
            "._suffix.SuffixValidator",
            "._prefixsrc.PrefixsrcValidator",
            "._prefix.PrefixValidator",
            "._line.LineValidator",
            "._height.HeightValidator",
            "._formatsrc.FormatsrcValidator",
            "._format.FormatValidator",
            "._font.FontValidator",
            "._fill.FillValidator",
            "._alignsrc.AlignsrcValidator",
            "._align.AlignValidator",
        ],
    )
