import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yside import YsideValidator
    from ._ygap import YgapValidator
    from ._yaxes import YaxesValidator
    from ._xside import XsideValidator
    from ._xgap import XgapValidator
    from ._xaxes import XaxesValidator
    from ._subplots import SubplotsValidator
    from ._rows import RowsValidator
    from ._roworder import RoworderValidator
    from ._pattern import PatternValidator
    from ._domain import DomainValidator
    from ._columns import ColumnsValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yside.YsideValidator",
            "._ygap.YgapValidator",
            "._yaxes.YaxesValidator",
            "._xside.XsideValidator",
            "._xgap.XgapValidator",
            "._xaxes.XaxesValidator",
            "._subplots.SubplotsValidator",
            "._rows.RowsValidator",
            "._roworder.RoworderValidator",
            "._pattern.PatternValidator",
            "._domain.DomainValidator",
            "._columns.ColumnsValidator",
        ],
    )
