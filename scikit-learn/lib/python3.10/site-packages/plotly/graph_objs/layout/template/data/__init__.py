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
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
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
