import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yshift import YshiftValidator
    from ._yref import YrefValidator
    from ._yclick import YclickValidator
    from ._yanchor import YanchorValidator
    from ._y import YValidator
    from ._xshift import XshiftValidator
    from ._xref import XrefValidator
    from ._xclick import XclickValidator
    from ._xanchor import XanchorValidator
    from ._x import XValidator
    from ._width import WidthValidator
    from ._visible import VisibleValidator
    from ._valign import ValignValidator
    from ._textangle import TextangleValidator
    from ._text import TextValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._startstandoff import StartstandoffValidator
    from ._startarrowsize import StartarrowsizeValidator
    from ._startarrowhead import StartarrowheadValidator
    from ._standoff import StandoffValidator
    from ._showarrow import ShowarrowValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._hovertext import HovertextValidator
    from ._hoverlabel import HoverlabelValidator
    from ._height import HeightValidator
    from ._font import FontValidator
    from ._clicktoshow import ClicktoshowValidator
    from ._captureevents import CaptureeventsValidator
    from ._borderwidth import BorderwidthValidator
    from ._borderpad import BorderpadValidator
    from ._bordercolor import BordercolorValidator
    from ._bgcolor import BgcolorValidator
    from ._ayref import AyrefValidator
    from ._ay import AyValidator
    from ._axref import AxrefValidator
    from ._ax import AxValidator
    from ._arrowwidth import ArrowwidthValidator
    from ._arrowsize import ArrowsizeValidator
    from ._arrowside import ArrowsideValidator
    from ._arrowhead import ArrowheadValidator
    from ._arrowcolor import ArrowcolorValidator
    from ._align import AlignValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yshift.YshiftValidator",
            "._yref.YrefValidator",
            "._yclick.YclickValidator",
            "._yanchor.YanchorValidator",
            "._y.YValidator",
            "._xshift.XshiftValidator",
            "._xref.XrefValidator",
            "._xclick.XclickValidator",
            "._xanchor.XanchorValidator",
            "._x.XValidator",
            "._width.WidthValidator",
            "._visible.VisibleValidator",
            "._valign.ValignValidator",
            "._textangle.TextangleValidator",
            "._text.TextValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._startstandoff.StartstandoffValidator",
            "._startarrowsize.StartarrowsizeValidator",
            "._startarrowhead.StartarrowheadValidator",
            "._standoff.StandoffValidator",
            "._showarrow.ShowarrowValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._hovertext.HovertextValidator",
            "._hoverlabel.HoverlabelValidator",
            "._height.HeightValidator",
            "._font.FontValidator",
            "._clicktoshow.ClicktoshowValidator",
            "._captureevents.CaptureeventsValidator",
            "._borderwidth.BorderwidthValidator",
            "._borderpad.BorderpadValidator",
            "._bordercolor.BordercolorValidator",
            "._bgcolor.BgcolorValidator",
            "._ayref.AyrefValidator",
            "._ay.AyValidator",
            "._axref.AxrefValidator",
            "._ax.AxValidator",
            "._arrowwidth.ArrowwidthValidator",
            "._arrowsize.ArrowsizeValidator",
            "._arrowside.ArrowsideValidator",
            "._arrowhead.ArrowheadValidator",
            "._arrowcolor.ArrowcolorValidator",
            "._align.AlignValidator",
        ],
    )
