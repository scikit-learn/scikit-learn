import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._valuesrc import ValuesrcValidator
    from ._value import ValueValidator
    from ._targetsrc import TargetsrcValidator
    from ._target import TargetValidator
    from ._sourcesrc import SourcesrcValidator
    from ._source import SourceValidator
    from ._line import LineValidator
    from ._labelsrc import LabelsrcValidator
    from ._label import LabelValidator
    from ._hovertemplatesrc import HovertemplatesrcValidator
    from ._hovertemplate import HovertemplateValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverinfo import HoverinfoValidator
    from ._hovercolorsrc import HovercolorsrcValidator
    from ._hovercolor import HovercolorValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._colorsrc import ColorsrcValidator
    from ._colorscaledefaults import ColorscaledefaultsValidator
    from ._colorscales import ColorscalesValidator
    from ._color import ColorValidator
    from ._arrowlen import ArrowlenValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._valuesrc.ValuesrcValidator",
            "._value.ValueValidator",
            "._targetsrc.TargetsrcValidator",
            "._target.TargetValidator",
            "._sourcesrc.SourcesrcValidator",
            "._source.SourceValidator",
            "._line.LineValidator",
            "._labelsrc.LabelsrcValidator",
            "._label.LabelValidator",
            "._hovertemplatesrc.HovertemplatesrcValidator",
            "._hovertemplate.HovertemplateValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverinfo.HoverinfoValidator",
            "._hovercolorsrc.HovercolorsrcValidator",
            "._hovercolor.HovercolorValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._colorsrc.ColorsrcValidator",
            "._colorscaledefaults.ColorscaledefaultsValidator",
            "._colorscales.ColorscalesValidator",
            "._color.ColorValidator",
            "._arrowlen.ArrowlenValidator",
        ],
    )
