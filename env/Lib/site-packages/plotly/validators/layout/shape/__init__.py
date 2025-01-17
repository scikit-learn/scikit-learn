import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._ysizemode import YsizemodeValidator
    from ._yref import YrefValidator
    from ._yanchor import YanchorValidator
    from ._y1shift import Y1ShiftValidator
    from ._y1 import Y1Validator
    from ._y0shift import Y0ShiftValidator
    from ._y0 import Y0Validator
    from ._xsizemode import XsizemodeValidator
    from ._xref import XrefValidator
    from ._xanchor import XanchorValidator
    from ._x1shift import X1ShiftValidator
    from ._x1 import X1Validator
    from ._x0shift import X0ShiftValidator
    from ._x0 import X0Validator
    from ._visible import VisibleValidator
    from ._type import TypeValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._showlegend import ShowlegendValidator
    from ._path import PathValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._line import LineValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legendgroup import LegendgroupValidator
    from ._legend import LegendValidator
    from ._layer import LayerValidator
    from ._label import LabelValidator
    from ._fillrule import FillruleValidator
    from ._fillcolor import FillcolorValidator
    from ._editable import EditableValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._ysizemode.YsizemodeValidator",
            "._yref.YrefValidator",
            "._yanchor.YanchorValidator",
            "._y1shift.Y1ShiftValidator",
            "._y1.Y1Validator",
            "._y0shift.Y0ShiftValidator",
            "._y0.Y0Validator",
            "._xsizemode.XsizemodeValidator",
            "._xref.XrefValidator",
            "._xanchor.XanchorValidator",
            "._x1shift.X1ShiftValidator",
            "._x1.X1Validator",
            "._x0shift.X0ShiftValidator",
            "._x0.X0Validator",
            "._visible.VisibleValidator",
            "._type.TypeValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._showlegend.ShowlegendValidator",
            "._path.PathValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._line.LineValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legendgroup.LegendgroupValidator",
            "._legend.LegendValidator",
            "._layer.LayerValidator",
            "._label.LabelValidator",
            "._fillrule.FillruleValidator",
            "._fillcolor.FillcolorValidator",
            "._editable.EditableValidator",
        ],
    )
