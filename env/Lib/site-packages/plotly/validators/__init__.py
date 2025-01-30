import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._waterfall import WaterfallValidator
    from ._volume import VolumeValidator
    from ._violin import ViolinValidator
    from ._treemap import TreemapValidator
    from ._table import TableValidator
    from ._surface import SurfaceValidator
    from ._sunburst import SunburstValidator
    from ._streamtube import StreamtubeValidator
    from ._splom import SplomValidator
    from ._scatterternary import ScatterternaryValidator
    from ._scattersmith import ScattersmithValidator
    from ._scatterpolargl import ScatterpolarglValidator
    from ._scatterpolar import ScatterpolarValidator
    from ._scattermapbox import ScattermapboxValidator
    from ._scattermap import ScattermapValidator
    from ._scattergl import ScatterglValidator
    from ._scattergeo import ScattergeoValidator
    from ._scattercarpet import ScattercarpetValidator
    from ._scatter3d import Scatter3DValidator
    from ._scatter import ScatterValidator
    from ._sankey import SankeyValidator
    from ._pointcloud import PointcloudValidator
    from ._pie import PieValidator
    from ._parcoords import ParcoordsValidator
    from ._parcats import ParcatsValidator
    from ._ohlc import OhlcValidator
    from ._mesh3d import Mesh3DValidator
    from ._isosurface import IsosurfaceValidator
    from ._indicator import IndicatorValidator
    from ._image import ImageValidator
    from ._icicle import IcicleValidator
    from ._histogram2dcontour import Histogram2DcontourValidator
    from ._histogram2d import Histogram2DValidator
    from ._histogram import HistogramValidator
    from ._heatmapgl import HeatmapglValidator
    from ._heatmap import HeatmapValidator
    from ._funnelarea import FunnelareaValidator
    from ._funnel import FunnelValidator
    from ._densitymapbox import DensitymapboxValidator
    from ._densitymap import DensitymapValidator
    from ._contourcarpet import ContourcarpetValidator
    from ._contour import ContourValidator
    from ._cone import ConeValidator
    from ._choroplethmapbox import ChoroplethmapboxValidator
    from ._choroplethmap import ChoroplethmapValidator
    from ._choropleth import ChoroplethValidator
    from ._carpet import CarpetValidator
    from ._candlestick import CandlestickValidator
    from ._box import BoxValidator
    from ._barpolar import BarpolarValidator
    from ._bar import BarValidator
    from ._layout import LayoutValidator
    from ._frames import FramesValidator
    from ._data import DataValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._waterfall.WaterfallValidator",
            "._volume.VolumeValidator",
            "._violin.ViolinValidator",
            "._treemap.TreemapValidator",
            "._table.TableValidator",
            "._surface.SurfaceValidator",
            "._sunburst.SunburstValidator",
            "._streamtube.StreamtubeValidator",
            "._splom.SplomValidator",
            "._scatterternary.ScatterternaryValidator",
            "._scattersmith.ScattersmithValidator",
            "._scatterpolargl.ScatterpolarglValidator",
            "._scatterpolar.ScatterpolarValidator",
            "._scattermapbox.ScattermapboxValidator",
            "._scattermap.ScattermapValidator",
            "._scattergl.ScatterglValidator",
            "._scattergeo.ScattergeoValidator",
            "._scattercarpet.ScattercarpetValidator",
            "._scatter3d.Scatter3DValidator",
            "._scatter.ScatterValidator",
            "._sankey.SankeyValidator",
            "._pointcloud.PointcloudValidator",
            "._pie.PieValidator",
            "._parcoords.ParcoordsValidator",
            "._parcats.ParcatsValidator",
            "._ohlc.OhlcValidator",
            "._mesh3d.Mesh3DValidator",
            "._isosurface.IsosurfaceValidator",
            "._indicator.IndicatorValidator",
            "._image.ImageValidator",
            "._icicle.IcicleValidator",
            "._histogram2dcontour.Histogram2DcontourValidator",
            "._histogram2d.Histogram2DValidator",
            "._histogram.HistogramValidator",
            "._heatmapgl.HeatmapglValidator",
            "._heatmap.HeatmapValidator",
            "._funnelarea.FunnelareaValidator",
            "._funnel.FunnelValidator",
            "._densitymapbox.DensitymapboxValidator",
            "._densitymap.DensitymapValidator",
            "._contourcarpet.ContourcarpetValidator",
            "._contour.ContourValidator",
            "._cone.ConeValidator",
            "._choroplethmapbox.ChoroplethmapboxValidator",
            "._choroplethmap.ChoroplethmapValidator",
            "._choropleth.ChoroplethValidator",
            "._carpet.CarpetValidator",
            "._candlestick.CandlestickValidator",
            "._box.BoxValidator",
            "._barpolar.BarpolarValidator",
            "._bar.BarValidator",
            "._layout.LayoutValidator",
            "._frames.FramesValidator",
            "._data.DataValidator",
        ],
    )
