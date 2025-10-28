import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._bar import Bar
    from ._barpolar import Barpolar
    from ._box import Box
    from ._candlestick import Candlestick
    from ._carpet import Carpet
    from ._choropleth import Choropleth
    from ._choroplethmap import Choroplethmap
    from ._choroplethmapbox import Choroplethmapbox
    from ._cone import Cone
    from ._contour import Contour
    from ._contourcarpet import Contourcarpet
    from ._densitymap import Densitymap
    from ._densitymapbox import Densitymapbox
    from ._deprecations import AngularAxis
    from ._deprecations import Annotation
    from ._deprecations import Annotations
    from ._deprecations import ColorBar
    from ._deprecations import Contours
    from ._deprecations import Data
    from ._deprecations import ErrorX
    from ._deprecations import ErrorY
    from ._deprecations import ErrorZ
    from ._deprecations import Font
    from ._deprecations import Frames
    from ._deprecations import Histogram2dcontour
    from ._deprecations import Legend
    from ._deprecations import Line
    from ._deprecations import Margin
    from ._deprecations import Marker
    from ._deprecations import RadialAxis
    from ._deprecations import Scene
    from ._deprecations import Stream
    from ._deprecations import Trace
    from ._deprecations import XAxis
    from ._deprecations import XBins
    from ._deprecations import YAxis
    from ._deprecations import YBins
    from ._deprecations import ZAxis
    from ._figure import Figure
    from ._frame import Frame
    from ._funnel import Funnel
    from ._funnelarea import Funnelarea
    from ._heatmap import Heatmap
    from ._histogram import Histogram
    from ._histogram2d import Histogram2d
    from ._histogram2dcontour import Histogram2dContour
    from ._icicle import Icicle
    from ._image import Image
    from ._indicator import Indicator
    from ._isosurface import Isosurface
    from ._layout import Layout
    from ._mesh3d import Mesh3d
    from ._ohlc import Ohlc
    from ._parcats import Parcats
    from ._parcoords import Parcoords
    from ._pie import Pie
    from ._sankey import Sankey
    from ._scatter import Scatter
    from ._scatter3d import Scatter3d
    from ._scattercarpet import Scattercarpet
    from ._scattergeo import Scattergeo
    from ._scattergl import Scattergl
    from ._scattermap import Scattermap
    from ._scattermapbox import Scattermapbox
    from ._scatterpolar import Scatterpolar
    from ._scatterpolargl import Scatterpolargl
    from ._scattersmith import Scattersmith
    from ._scatterternary import Scatterternary
    from ._splom import Splom
    from ._streamtube import Streamtube
    from ._sunburst import Sunburst
    from ._surface import Surface
    from ._table import Table
    from ._treemap import Treemap
    from ._violin import Violin
    from ._volume import Volume
    from ._waterfall import Waterfall
    from . import bar
    from . import barpolar
    from . import box
    from . import candlestick
    from . import carpet
    from . import choropleth
    from . import choroplethmap
    from . import choroplethmapbox
    from . import cone
    from . import contour
    from . import contourcarpet
    from . import densitymap
    from . import densitymapbox
    from . import funnel
    from . import funnelarea
    from . import heatmap
    from . import histogram
    from . import histogram2d
    from . import histogram2dcontour
    from . import icicle
    from . import image
    from . import indicator
    from . import isosurface
    from . import layout
    from . import mesh3d
    from . import ohlc
    from . import parcats
    from . import parcoords
    from . import pie
    from . import sankey
    from . import scatter
    from . import scatter3d
    from . import scattercarpet
    from . import scattergeo
    from . import scattergl
    from . import scattermap
    from . import scattermapbox
    from . import scatterpolar
    from . import scatterpolargl
    from . import scattersmith
    from . import scatterternary
    from . import splom
    from . import streamtube
    from . import sunburst
    from . import surface
    from . import table
    from . import treemap
    from . import violin
    from . import volume
    from . import waterfall
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [
            ".bar",
            ".barpolar",
            ".box",
            ".candlestick",
            ".carpet",
            ".choropleth",
            ".choroplethmap",
            ".choroplethmapbox",
            ".cone",
            ".contour",
            ".contourcarpet",
            ".densitymap",
            ".densitymapbox",
            ".funnel",
            ".funnelarea",
            ".heatmap",
            ".histogram",
            ".histogram2d",
            ".histogram2dcontour",
            ".icicle",
            ".image",
            ".indicator",
            ".isosurface",
            ".layout",
            ".mesh3d",
            ".ohlc",
            ".parcats",
            ".parcoords",
            ".pie",
            ".sankey",
            ".scatter",
            ".scatter3d",
            ".scattercarpet",
            ".scattergeo",
            ".scattergl",
            ".scattermap",
            ".scattermapbox",
            ".scatterpolar",
            ".scatterpolargl",
            ".scattersmith",
            ".scatterternary",
            ".splom",
            ".streamtube",
            ".sunburst",
            ".surface",
            ".table",
            ".treemap",
            ".violin",
            ".volume",
            ".waterfall",
        ],
        [
            "._bar.Bar",
            "._barpolar.Barpolar",
            "._box.Box",
            "._candlestick.Candlestick",
            "._carpet.Carpet",
            "._choropleth.Choropleth",
            "._choroplethmap.Choroplethmap",
            "._choroplethmapbox.Choroplethmapbox",
            "._cone.Cone",
            "._contour.Contour",
            "._contourcarpet.Contourcarpet",
            "._densitymap.Densitymap",
            "._densitymapbox.Densitymapbox",
            "._deprecations.AngularAxis",
            "._deprecations.Annotation",
            "._deprecations.Annotations",
            "._deprecations.ColorBar",
            "._deprecations.Contours",
            "._deprecations.Data",
            "._deprecations.ErrorX",
            "._deprecations.ErrorY",
            "._deprecations.ErrorZ",
            "._deprecations.Font",
            "._deprecations.Frames",
            "._deprecations.Histogram2dcontour",
            "._deprecations.Legend",
            "._deprecations.Line",
            "._deprecations.Margin",
            "._deprecations.Marker",
            "._deprecations.RadialAxis",
            "._deprecations.Scene",
            "._deprecations.Stream",
            "._deprecations.Trace",
            "._deprecations.XAxis",
            "._deprecations.XBins",
            "._deprecations.YAxis",
            "._deprecations.YBins",
            "._deprecations.ZAxis",
            "._figure.Figure",
            "._frame.Frame",
            "._funnel.Funnel",
            "._funnelarea.Funnelarea",
            "._heatmap.Heatmap",
            "._histogram.Histogram",
            "._histogram2d.Histogram2d",
            "._histogram2dcontour.Histogram2dContour",
            "._icicle.Icicle",
            "._image.Image",
            "._indicator.Indicator",
            "._isosurface.Isosurface",
            "._layout.Layout",
            "._mesh3d.Mesh3d",
            "._ohlc.Ohlc",
            "._parcats.Parcats",
            "._parcoords.Parcoords",
            "._pie.Pie",
            "._sankey.Sankey",
            "._scatter.Scatter",
            "._scatter3d.Scatter3d",
            "._scattercarpet.Scattercarpet",
            "._scattergeo.Scattergeo",
            "._scattergl.Scattergl",
            "._scattermap.Scattermap",
            "._scattermapbox.Scattermapbox",
            "._scatterpolar.Scatterpolar",
            "._scatterpolargl.Scatterpolargl",
            "._scattersmith.Scattersmith",
            "._scatterternary.Scatterternary",
            "._splom.Splom",
            "._streamtube.Streamtube",
            "._sunburst.Sunburst",
            "._surface.Surface",
            "._table.Table",
            "._treemap.Treemap",
            "._violin.Violin",
            "._volume.Volume",
            "._waterfall.Waterfall",
        ],
    )


if sys.version_info < (3, 7) or TYPE_CHECKING:
    try:
        import ipywidgets as _ipywidgets
        from packaging.version import Version as _Version

        if _Version(_ipywidgets.__version__) >= _Version("7.0.0"):
            from ..graph_objs._figurewidget import FigureWidget
        else:
            raise ImportError()
    except Exception:
        from ..missing_anywidget import FigureWidget
else:
    __all__.append("FigureWidget")
    orig_getattr = __getattr__

    def __getattr__(import_name):
        if import_name == "FigureWidget":
            try:
                import ipywidgets
                from packaging.version import Version

                if Version(ipywidgets.__version__) >= Version("7.0.0"):
                    from ..graph_objs._figurewidget import FigureWidget

                    return FigureWidget
                else:
                    raise ImportError()
            except Exception:
                from ..missing_anywidget import FigureWidget

                return FigureWidget
            else:
                raise ImportError()

        return orig_getattr(import_name)
