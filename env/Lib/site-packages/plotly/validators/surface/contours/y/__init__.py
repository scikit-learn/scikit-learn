import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._width import WidthValidator
    from ._usecolormap import UsecolormapValidator
    from ._start import StartValidator
    from ._size import SizeValidator
    from ._show import ShowValidator
    from ._project import ProjectValidator
    from ._highlightwidth import HighlightwidthValidator
    from ._highlightcolor import HighlightcolorValidator
    from ._highlight import HighlightValidator
    from ._end import EndValidator
    from ._color import ColorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._width.WidthValidator",
            "._usecolormap.UsecolormapValidator",
            "._start.StartValidator",
            "._size.SizeValidator",
            "._show.ShowValidator",
            "._project.ProjectValidator",
            "._highlightwidth.HighlightwidthValidator",
            "._highlightcolor.HighlightcolorValidator",
            "._highlight.HighlightValidator",
            "._end.EndValidator",
            "._color.ColorValidator",
        ],
    )
