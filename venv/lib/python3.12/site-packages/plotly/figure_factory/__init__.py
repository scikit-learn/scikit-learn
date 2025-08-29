# ruff: noqa: E402

from plotly import optional_imports

# Require that numpy exists for figure_factory
np = optional_imports.get_module("numpy")
if np is None:
    raise ImportError(
        """\
The figure factory module requires the numpy package"""
    )


from plotly.figure_factory._2d_density import create_2d_density
from plotly.figure_factory._annotated_heatmap import create_annotated_heatmap
from plotly.figure_factory._bullet import create_bullet
from plotly.figure_factory._candlestick import create_candlestick
from plotly.figure_factory._dendrogram import create_dendrogram
from plotly.figure_factory._distplot import create_distplot
from plotly.figure_factory._facet_grid import create_facet_grid
from plotly.figure_factory._gantt import create_gantt
from plotly.figure_factory._ohlc import create_ohlc
from plotly.figure_factory._quiver import create_quiver
from plotly.figure_factory._scatterplot import create_scatterplotmatrix
from plotly.figure_factory._streamline import create_streamline
from plotly.figure_factory._table import create_table
from plotly.figure_factory._trisurf import create_trisurf
from plotly.figure_factory._violin import create_violin

if optional_imports.get_module("pandas") is not None:
    from plotly.figure_factory._county_choropleth import create_choropleth
    from plotly.figure_factory._hexbin_mapbox import create_hexbin_mapbox
else:

    def create_choropleth(*args, **kwargs):
        raise ImportError("Please install pandas to use `create_choropleth`")

    def create_hexbin_mapbox(*args, **kwargs):
        raise ImportError("Please install pandas to use `create_hexbin_mapbox`")


if optional_imports.get_module("skimage") is not None:
    from plotly.figure_factory._ternary_contour import create_ternary_contour
else:

    def create_ternary_contour(*args, **kwargs):
        raise ImportError("Please install scikit-image to use `create_ternary_contour`")


__all__ = [
    "create_2d_density",
    "create_annotated_heatmap",
    "create_bullet",
    "create_candlestick",
    "create_choropleth",
    "create_dendrogram",
    "create_distplot",
    "create_facet_grid",
    "create_gantt",
    "create_hexbin_mapbox",
    "create_ohlc",
    "create_quiver",
    "create_scatterplotmatrix",
    "create_streamline",
    "create_table",
    "create_ternary_contour",
    "create_trisurf",
    "create_violin",
]
