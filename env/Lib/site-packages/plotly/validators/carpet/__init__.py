import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._zorder import ZorderValidator
    from ._ysrc import YsrcValidator
    from ._yaxis import YaxisValidator
    from ._y import YValidator
    from ._xsrc import XsrcValidator
    from ._xaxis import XaxisValidator
    from ._x import XValidator
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._uid import UidValidator
    from ._stream import StreamValidator
    from ._opacity import OpacityValidator
    from ._name import NameValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._legendwidth import LegendwidthValidator
    from ._legendrank import LegendrankValidator
    from ._legendgrouptitle import LegendgrouptitleValidator
    from ._legend import LegendValidator
    from ._idssrc import IdssrcValidator
    from ._ids import IdsValidator
    from ._font import FontValidator
    from ._db import DbValidator
    from ._da import DaValidator
    from ._customdatasrc import CustomdatasrcValidator
    from ._customdata import CustomdataValidator
    from ._color import ColorValidator
    from ._cheaterslope import CheaterslopeValidator
    from ._carpet import CarpetValidator
    from ._bsrc import BsrcValidator
    from ._baxis import BaxisValidator
    from ._b0 import B0Validator
    from ._b import BValidator
    from ._asrc import AsrcValidator
    from ._aaxis import AaxisValidator
    from ._a0 import A0Validator
    from ._a import AValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._zorder.ZorderValidator",
            "._ysrc.YsrcValidator",
            "._yaxis.YaxisValidator",
            "._y.YValidator",
            "._xsrc.XsrcValidator",
            "._xaxis.XaxisValidator",
            "._x.XValidator",
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._uid.UidValidator",
            "._stream.StreamValidator",
            "._opacity.OpacityValidator",
            "._name.NameValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._legendwidth.LegendwidthValidator",
            "._legendrank.LegendrankValidator",
            "._legendgrouptitle.LegendgrouptitleValidator",
            "._legend.LegendValidator",
            "._idssrc.IdssrcValidator",
            "._ids.IdsValidator",
            "._font.FontValidator",
            "._db.DbValidator",
            "._da.DaValidator",
            "._customdatasrc.CustomdatasrcValidator",
            "._customdata.CustomdataValidator",
            "._color.ColorValidator",
            "._cheaterslope.CheaterslopeValidator",
            "._carpet.CarpetValidator",
            "._bsrc.BsrcValidator",
            "._baxis.BaxisValidator",
            "._b0.B0Validator",
            "._b.BValidator",
            "._asrc.AsrcValidator",
            "._aaxis.AaxisValidator",
            "._a0.A0Validator",
            "._a.AValidator",
        ],
    )
