# ruff: noqa: F401
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..graph_objs import Waterfall
    from ..graph_objs import Volume
    from ..graph_objs import Violin
    from ..graph_objs import Treemap
    from ..graph_objs import Table
    from ..graph_objs import Surface
    from ..graph_objs import Sunburst
    from ..graph_objs import Streamtube
    from ..graph_objs import Splom
    from ..graph_objs import Scatterternary
    from ..graph_objs import Scattersmith
    from ..graph_objs import Scatterpolargl
    from ..graph_objs import Scatterpolar
    from ..graph_objs import Scattermapbox
    from ..graph_objs import Scattermap
    from ..graph_objs import Scattergl
    from ..graph_objs import Scattergeo
    from ..graph_objs import Scattercarpet
    from ..graph_objs import Scatter3d
    from ..graph_objs import Scatter
    from ..graph_objs import Sankey
    from ..graph_objs import Pie
    from ..graph_objs import Parcoords
    from ..graph_objs import Parcats
    from ..graph_objs import Ohlc
    from ..graph_objs import Mesh3d
    from ..graph_objs import Isosurface
    from ..graph_objs import Indicator
    from ..graph_objs import Image
    from ..graph_objs import Icicle
    from ..graph_objs import Histogram2dContour
    from ..graph_objs import Histogram2d
    from ..graph_objs import Histogram
    from ..graph_objs import Heatmap
    from ..graph_objs import Funnelarea
    from ..graph_objs import Funnel
    from ..graph_objs import Densitymapbox
    from ..graph_objs import Densitymap
    from ..graph_objs import Contourcarpet
    from ..graph_objs import Contour
    from ..graph_objs import Cone
    from ..graph_objs import Choroplethmapbox
    from ..graph_objs import Choroplethmap
    from ..graph_objs import Choropleth
    from ..graph_objs import Carpet
    from ..graph_objs import Candlestick
    from ..graph_objs import Box
    from ..graph_objs import Barpolar
    from ..graph_objs import Bar
    from ..graph_objs import Layout
    from ..graph_objs import Frame
    from ..graph_objs import Figure
    from ..graph_objs import Data
    from ..graph_objs import Annotations
    from ..graph_objs import Frames
    from ..graph_objs import AngularAxis
    from ..graph_objs import Annotation
    from ..graph_objs import ColorBar
    from ..graph_objs import Contours
    from ..graph_objs import ErrorX
    from ..graph_objs import ErrorY
    from ..graph_objs import ErrorZ
    from ..graph_objs import Font
    from ..graph_objs import Legend
    from ..graph_objs import Line
    from ..graph_objs import Margin
    from ..graph_objs import Marker
    from ..graph_objs import RadialAxis
    from ..graph_objs import Scene
    from ..graph_objs import Stream
    from ..graph_objs import XAxis
    from ..graph_objs import YAxis
    from ..graph_objs import ZAxis
    from ..graph_objs import XBins
    from ..graph_objs import YBins
    from ..graph_objs import Trace
    from ..graph_objs import Histogram2dcontour
    from ..graph_objs import waterfall
    from ..graph_objs import volume
    from ..graph_objs import violin
    from ..graph_objs import treemap
    from ..graph_objs import table
    from ..graph_objs import surface
    from ..graph_objs import sunburst
    from ..graph_objs import streamtube
    from ..graph_objs import splom
    from ..graph_objs import scatterternary
    from ..graph_objs import scattersmith
    from ..graph_objs import scatterpolargl
    from ..graph_objs import scatterpolar
    from ..graph_objs import scattermapbox
    from ..graph_objs import scattermap
    from ..graph_objs import scattergl
    from ..graph_objs import scattergeo
    from ..graph_objs import scattercarpet
    from ..graph_objs import scatter3d
    from ..graph_objs import scatter
    from ..graph_objs import sankey
    from ..graph_objs import pie
    from ..graph_objs import parcoords
    from ..graph_objs import parcats
    from ..graph_objs import ohlc
    from ..graph_objs import mesh3d
    from ..graph_objs import isosurface
    from ..graph_objs import indicator
    from ..graph_objs import image
    from ..graph_objs import icicle
    from ..graph_objs import histogram2dcontour
    from ..graph_objs import histogram2d
    from ..graph_objs import histogram
    from ..graph_objs import heatmap
    from ..graph_objs import funnelarea
    from ..graph_objs import funnel
    from ..graph_objs import densitymapbox
    from ..graph_objs import densitymap
    from ..graph_objs import contourcarpet
    from ..graph_objs import contour
    from ..graph_objs import cone
    from ..graph_objs import choroplethmapbox
    from ..graph_objs import choroplethmap
    from ..graph_objs import choropleth
    from ..graph_objs import carpet
    from ..graph_objs import candlestick
    from ..graph_objs import box
    from ..graph_objs import barpolar
    from ..graph_objs import bar
    from ..graph_objs import layout
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [
            "..graph_objs.waterfall",
            "..graph_objs.volume",
            "..graph_objs.violin",
            "..graph_objs.treemap",
            "..graph_objs.table",
            "..graph_objs.surface",
            "..graph_objs.sunburst",
            "..graph_objs.streamtube",
            "..graph_objs.splom",
            "..graph_objs.scatterternary",
            "..graph_objs.scattersmith",
            "..graph_objs.scatterpolargl",
            "..graph_objs.scatterpolar",
            "..graph_objs.scattermapbox",
            "..graph_objs.scattermap",
            "..graph_objs.scattergl",
            "..graph_objs.scattergeo",
            "..graph_objs.scattercarpet",
            "..graph_objs.scatter3d",
            "..graph_objs.scatter",
            "..graph_objs.sankey",
            "..graph_objs.pie",
            "..graph_objs.parcoords",
            "..graph_objs.parcats",
            "..graph_objs.ohlc",
            "..graph_objs.mesh3d",
            "..graph_objs.isosurface",
            "..graph_objs.indicator",
            "..graph_objs.image",
            "..graph_objs.icicle",
            "..graph_objs.histogram2dcontour",
            "..graph_objs.histogram2d",
            "..graph_objs.histogram",
            "..graph_objs.heatmap",
            "..graph_objs.funnelarea",
            "..graph_objs.funnel",
            "..graph_objs.densitymapbox",
            "..graph_objs.densitymap",
            "..graph_objs.contourcarpet",
            "..graph_objs.contour",
            "..graph_objs.cone",
            "..graph_objs.choroplethmapbox",
            "..graph_objs.choroplethmap",
            "..graph_objs.choropleth",
            "..graph_objs.carpet",
            "..graph_objs.candlestick",
            "..graph_objs.box",
            "..graph_objs.barpolar",
            "..graph_objs.bar",
            "..graph_objs.layout",
        ],
        [
            "..graph_objs.Waterfall",
            "..graph_objs.Volume",
            "..graph_objs.Violin",
            "..graph_objs.Treemap",
            "..graph_objs.Table",
            "..graph_objs.Surface",
            "..graph_objs.Sunburst",
            "..graph_objs.Streamtube",
            "..graph_objs.Splom",
            "..graph_objs.Scatterternary",
            "..graph_objs.Scattersmith",
            "..graph_objs.Scatterpolargl",
            "..graph_objs.Scatterpolar",
            "..graph_objs.Scattermapbox",
            "..graph_objs.Scattermap",
            "..graph_objs.Scattergl",
            "..graph_objs.Scattergeo",
            "..graph_objs.Scattercarpet",
            "..graph_objs.Scatter3d",
            "..graph_objs.Scatter",
            "..graph_objs.Sankey",
            "..graph_objs.Pie",
            "..graph_objs.Parcoords",
            "..graph_objs.Parcats",
            "..graph_objs.Ohlc",
            "..graph_objs.Mesh3d",
            "..graph_objs.Isosurface",
            "..graph_objs.Indicator",
            "..graph_objs.Image",
            "..graph_objs.Icicle",
            "..graph_objs.Histogram2dContour",
            "..graph_objs.Histogram2d",
            "..graph_objs.Histogram",
            "..graph_objs.Heatmap",
            "..graph_objs.Funnelarea",
            "..graph_objs.Funnel",
            "..graph_objs.Densitymapbox",
            "..graph_objs.Densitymap",
            "..graph_objs.Contourcarpet",
            "..graph_objs.Contour",
            "..graph_objs.Cone",
            "..graph_objs.Choroplethmapbox",
            "..graph_objs.Choroplethmap",
            "..graph_objs.Choropleth",
            "..graph_objs.Carpet",
            "..graph_objs.Candlestick",
            "..graph_objs.Box",
            "..graph_objs.Barpolar",
            "..graph_objs.Bar",
            "..graph_objs.Layout",
            "..graph_objs.Frame",
            "..graph_objs.Figure",
            "..graph_objs.Data",
            "..graph_objs.Annotations",
            "..graph_objs.Frames",
            "..graph_objs.AngularAxis",
            "..graph_objs.Annotation",
            "..graph_objs.ColorBar",
            "..graph_objs.Contours",
            "..graph_objs.ErrorX",
            "..graph_objs.ErrorY",
            "..graph_objs.ErrorZ",
            "..graph_objs.Font",
            "..graph_objs.Legend",
            "..graph_objs.Line",
            "..graph_objs.Margin",
            "..graph_objs.Marker",
            "..graph_objs.RadialAxis",
            "..graph_objs.Scene",
            "..graph_objs.Stream",
            "..graph_objs.XAxis",
            "..graph_objs.YAxis",
            "..graph_objs.ZAxis",
            "..graph_objs.XBins",
            "..graph_objs.YBins",
            "..graph_objs.Trace",
            "..graph_objs.Histogram2dcontour",
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
