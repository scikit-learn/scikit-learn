import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._value import ValueValidator
    from ._type import TypeValidator
    from ._start import StartValidator
    from ._size import SizeValidator
    from ._showlines import ShowlinesValidator
    from ._showlabels import ShowlabelsValidator
    from ._operation import OperationValidator
    from ._labelformat import LabelformatValidator
    from ._labelfont import LabelfontValidator
    from ._end import EndValidator
    from ._coloring import ColoringValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._value.ValueValidator",
            "._type.TypeValidator",
            "._start.StartValidator",
            "._size.SizeValidator",
            "._showlines.ShowlinesValidator",
            "._showlabels.ShowlabelsValidator",
            "._operation.OperationValidator",
            "._labelformat.LabelformatValidator",
            "._labelfont.LabelfontValidator",
            "._end.EndValidator",
            "._coloring.ColoringValidator",
        ],
    )
