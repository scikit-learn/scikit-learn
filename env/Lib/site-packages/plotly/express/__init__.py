"""
`plotly.express` is a terse, consistent, high-level wrapper around `plotly.graph_objects`
for rapid data exploration and figure generation. Learn more at https://plotly.com/python/plotly-express/
"""
from plotly import optional_imports

pd = optional_imports.get_module("pandas")
if pd is None:
    raise ImportError(
        """\
Plotly express requires pandas to be installed."""
    )

from ._imshow import imshow
from ._chart_types import (  # noqa: F401
    scatter,
    scatter_3d,
    scatter_polar,
    scatter_ternary,
    scatter_map,
    scatter_mapbox,
    scatter_geo,
    line,
    line_3d,
    line_polar,
    line_ternary,
    line_map,
    line_mapbox,
    line_geo,
    area,
    bar,
    timeline,
    bar_polar,
    violin,
    box,
    strip,
    histogram,
    ecdf,
    scatter_matrix,
    parallel_coordinates,
    parallel_categories,
    choropleth,
    density_contour,
    density_heatmap,
    pie,
    sunburst,
    treemap,
    icicle,
    funnel,
    funnel_area,
    choropleth_map,
    choropleth_mapbox,
    density_map,
    density_mapbox,
)


from ._core import (  # noqa: F401
    set_mapbox_access_token,
    defaults,
    get_trendline_results,
    NO_COLOR,
)

from ._special_inputs import IdentityMap, Constant, Range  # noqa: F401

from . import data, colors, trendline_functions  # noqa: F401

__all__ = [
    "scatter",
    "scatter_3d",
    "scatter_polar",
    "scatter_ternary",
    "scatter_map",
    "scatter_mapbox",
    "scatter_geo",
    "scatter_matrix",
    "density_contour",
    "density_heatmap",
    "density_map",
    "density_mapbox",
    "line",
    "line_3d",
    "line_polar",
    "line_ternary",
    "line_map",
    "line_mapbox",
    "line_geo",
    "parallel_coordinates",
    "parallel_categories",
    "area",
    "bar",
    "timeline",
    "bar_polar",
    "violin",
    "box",
    "strip",
    "histogram",
    "ecdf",
    "choropleth",
    "choropleth_map",
    "choropleth_mapbox",
    "pie",
    "sunburst",
    "treemap",
    "icicle",
    "funnel",
    "funnel_area",
    "imshow",
    "data",
    "colors",
    "trendline_functions",
    "set_mapbox_access_token",
    "get_trendline_results",
    "IdentityMap",
    "Constant",
    "Range",
    "NO_COLOR",
]
