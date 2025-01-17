import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._ysrc import YsrcValidator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._x import XValidator
    from ._thickness import ThicknessValidator
    from ._pad import PadValidator
    from ._line import LineValidator
    from ._labelsrc import LabelsrcValidator
    from ._label import LabelValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfo import HoverinfoValidator
    from ._groups import GroupsValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._colorsrc import ColorsrcValidator
    from ._color import ColorValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._ysrc.YsrcValidator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._x.XValidator",
            "._thickness.ThicknessValidator",
            "._pad.PadValidator",
            "._line.LineValidator",
            "._labelsrc.LabelsrcValidator",
            "._label.LabelValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfo.HoverinfoValidator",
            "._groups.GroupsValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._colorsrc.ColorsrcValidator",
            "._color.ColorValidator",
            "._align.AlignValidator",
        ],
    )
