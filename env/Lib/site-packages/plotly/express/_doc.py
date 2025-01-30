import inspect
from textwrap import TextWrapper

try:
    getfullargspec = inspect.getfullargspec
except AttributeError:  # python 2
    getfullargspec = inspect.getargspec


colref_type = "str or int or Series or array-like"
colref_desc = "Either a name of a column in `data_frame`, or a pandas Series or array_like object."
colref_list_type = "list of str or int, or Series or array-like"
colref_list_desc = (
    "Either names of columns in `data_frame`, or pandas Series, or array_like objects"
)

docs = dict(
    data_frame=[
        "DataFrame or array-like or dict",
        "This argument needs to be passed for column names (and not keyword names) to be used.",
        "Array-like and dict are transformed internally to a pandas DataFrame.",
        "Optional: if missing, a DataFrame gets constructed under the hood using the other arguments.",
    ],
    x=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the x axis in cartesian coordinates.",
    ],
    y=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the y axis in cartesian coordinates.",
    ],
    z=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the z axis in cartesian coordinates.",
    ],
    x_start=[
        colref_type,
        colref_desc,
        "(required)",
        "Values from this column or array_like are used to position marks along the x axis in cartesian coordinates.",
    ],
    x_end=[
        colref_type,
        colref_desc,
        "(required)",
        "Values from this column or array_like are used to position marks along the x axis in cartesian coordinates.",
    ],
    a=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the a axis in ternary coordinates.",
    ],
    b=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the b axis in ternary coordinates.",
    ],
    c=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the c axis in ternary coordinates.",
    ],
    r=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the radial axis in polar coordinates.",
    ],
    theta=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks along the angular axis in polar coordinates.",
    ],
    values=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to set values associated to sectors.",
    ],
    parents=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used as parents in sunburst and treemap charts.",
    ],
    ids=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to set ids of sectors",
    ],
    path=[
        colref_list_type,
        colref_list_desc,
        "List of columns names or columns of a rectangular dataframe defining the hierarchy of sectors, from root to leaves.",
        "An error is raised if path AND ids or parents is passed",
    ],
    lat=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks according to latitude on a map.",
    ],
    lon=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position marks according to longitude on a map.",
    ],
    locations=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are to be interpreted according to `locationmode` and mapped to longitude/latitude.",
    ],
    base=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to position the base of the bar.",
    ],
    dimensions=[
        colref_list_type,
        colref_list_desc,
        "Values from these columns are used for multidimensional visualization.",
    ],
    dimensions_max_cardinality=[
        "int (default 50)",
        "When `dimensions` is `None` and `data_frame` is provided, "
        "columns with more than this number of unique values are excluded from the output.",
        "Not used when `dimensions` is passed.",
    ],
    error_x=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size x-axis error bars.",
        "If `error_x_minus` is `None`, error bars will be symmetrical, otherwise `error_x` is used for the positive direction only.",
    ],
    error_x_minus=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size x-axis error bars in the negative direction.",
        "Ignored if `error_x` is `None`.",
    ],
    error_y=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size y-axis error bars.",
        "If `error_y_minus` is `None`, error bars will be symmetrical, otherwise `error_y` is used for the positive direction only.",
    ],
    error_y_minus=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size y-axis error bars in the negative direction.",
        "Ignored if `error_y` is `None`.",
    ],
    error_z=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size z-axis error bars.",
        "If `error_z_minus` is `None`, error bars will be symmetrical, otherwise `error_z` is used for the positive direction only.",
    ],
    error_z_minus=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to size z-axis error bars in the negative direction.",
        "Ignored if `error_z` is `None`.",
    ],
    color=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign color to marks.",
    ],
    opacity=["float", "Value between 0 and 1. Sets the opacity for markers."],
    line_dash=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign dash-patterns to lines.",
    ],
    line_group=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to group rows of `data_frame` into lines.",
    ],
    symbol=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign symbols to marks.",
    ],
    pattern_shape=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign pattern shapes to marks.",
    ],
    size=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign mark sizes.",
    ],
    radius=["int (default is 30)", "Sets the radius of influence of each point."],
    hover_name=[
        colref_type,
        colref_desc,
        "Values from this column or array_like appear in bold in the hover tooltip.",
    ],
    hover_data=[
        "str, or list of str or int, or Series or array-like, or dict",
        "Either a name or list of names of columns in `data_frame`, or pandas Series,",
        "or array_like objects",
        "or a dict with column names as keys, with values True (for default formatting)",
        "False (in order to remove this column from hover information),",
        "or a formatting string, for example ':.3f' or '|%a'",
        "or list-like data to appear in the hover tooltip",
        "or tuples with a bool or formatting string as first element,",
        "and list-like data to appear in hover as second element",
        "Values from these columns appear as extra data in the hover tooltip.",
    ],
    custom_data=[
        "str, or list of str or int, or Series or array-like",
        "Either name or list of names of columns in `data_frame`, or pandas Series, or array_like objects",
        "Values from these columns are extra data, to be used in widgets or Dash callbacks for example. This data is not user-visible but is included in events emitted by the figure (lasso selection etc.)",
    ],
    text=[
        colref_type,
        colref_desc,
        "Values from this column or array_like appear in the figure as text labels.",
    ],
    names=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used as labels for sectors.",
    ],
    locationmode=[
        "str",
        "One of 'ISO-3', 'USA-states', or 'country names'",
        "Determines the set of locations used to match entries in `locations` to regions on the map.",
    ],
    facet_row=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign marks to facetted subplots in the vertical direction.",
    ],
    facet_col=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign marks to facetted subplots in the horizontal direction.",
    ],
    facet_col_wrap=[
        "int",
        "Maximum number of facet columns.",
        "Wraps the column variable at this width, so that the column facets span multiple rows.",
        "Ignored if 0, and forced to 0 if `facet_row` or a `marginal` is set.",
    ],
    facet_row_spacing=[
        "float between 0 and 1",
        "Spacing between facet rows, in paper units. Default is 0.03 or 0.07 when facet_col_wrap is used.",
    ],
    facet_col_spacing=[
        "float between 0 and 1",
        "Spacing between facet columns, in paper units Default is 0.02.",
    ],
    animation_frame=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to assign marks to animation frames.",
    ],
    animation_group=[
        colref_type,
        colref_desc,
        "Values from this column or array_like are used to provide object-constancy across animation frames: rows with matching `animation_group`s will be treated as if they describe the same object in each frame.",
    ],
    symbol_sequence=[
        "list of str",
        "Strings should define valid plotly.js symbols.",
        "When `symbol` is set, values in that column are assigned symbols by cycling through `symbol_sequence` in the order described in `category_orders`, unless the value of `symbol` is a key in `symbol_map`.",
    ],
    symbol_map=[
        "dict with str keys and str values (default `{}`)",
        "String values should define plotly.js symbols",
        "Used to override `symbol_sequence` to assign a specific symbols to marks corresponding with specific values.",
        "Keys in `symbol_map` should be values in the column denoted by `symbol`.",
        "Alternatively, if the values of `symbol` are valid symbol names, the string `'identity'` may be passed to cause them to be used directly.",
    ],
    line_dash_map=[
        "dict with str keys and str values (default `{}`)",
        "Strings values define plotly.js dash-patterns.",
        "Used to override `line_dash_sequences` to assign a specific dash-patterns to lines corresponding with specific values.",
        "Keys in `line_dash_map` should be values in the column denoted by `line_dash`.",
        "Alternatively, if the values of `line_dash` are valid line-dash names, the string `'identity'` may be passed to cause them to be used directly.",
    ],
    line_dash_sequence=[
        "list of str",
        "Strings should define valid plotly.js dash-patterns.",
        "When `line_dash` is set, values in that column are assigned dash-patterns by cycling through `line_dash_sequence` in the order described in `category_orders`, unless the value of `line_dash` is a key in `line_dash_map`.",
    ],
    pattern_shape_map=[
        "dict with str keys and str values (default `{}`)",
        "Strings values define plotly.js patterns-shapes.",
        "Used to override `pattern_shape_sequences` to assign a specific patterns-shapes to lines corresponding with specific values.",
        "Keys in `pattern_shape_map` should be values in the column denoted by `pattern_shape`.",
        "Alternatively, if the values of `pattern_shape` are valid patterns-shapes names, the string `'identity'` may be passed to cause them to be used directly.",
    ],
    pattern_shape_sequence=[
        "list of str",
        "Strings should define valid plotly.js patterns-shapes.",
        "When `pattern_shape` is set, values in that column are assigned patterns-shapes by cycling through `pattern_shape_sequence` in the order described in `category_orders`, unless the value of `pattern_shape` is a key in `pattern_shape_map`.",
    ],
    color_discrete_sequence=[
        "list of str",
        "Strings should define valid CSS-colors.",
        "When `color` is set and the values in the corresponding column are not numeric, values in that column are assigned colors by cycling through `color_discrete_sequence` in the order described in `category_orders`, unless the value of `color` is a key in `color_discrete_map`.",
        "Various useful color sequences are available in the `plotly.express.colors` submodules, specifically `plotly.express.colors.qualitative`.",
    ],
    color_discrete_map=[
        "dict with str keys and str values (default `{}`)",
        "String values should define valid CSS-colors",
        "Used to override `color_discrete_sequence` to assign a specific colors to marks corresponding with specific values.",
        "Keys in `color_discrete_map` should be values in the column denoted by `color`.",
        "Alternatively, if the values of `color` are valid colors, the string `'identity'` may be passed to cause them to be used directly.",
    ],
    color_continuous_scale=[
        "list of str",
        "Strings should define valid CSS-colors",
        "This list is used to build a continuous color scale when the column denoted by `color` contains numeric data.",
        "Various useful color scales are available in the `plotly.express.colors` submodules, specifically `plotly.express.colors.sequential`, `plotly.express.colors.diverging` and `plotly.express.colors.cyclical`.",
    ],
    color_continuous_midpoint=[
        "number (default `None`)",
        "If set, computes the bounds of the continuous color scale to have the desired midpoint.",
        "Setting this value is recommended when using `plotly.express.colors.diverging` color scales as the inputs to `color_continuous_scale`.",
    ],
    size_max=["int (default `20`)", "Set the maximum mark size when using `size`."],
    markers=["boolean (default `False`)", "If `True`, markers are shown on lines."],
    lines=[
        "boolean (default `True`)",
        "If `False`, lines are not drawn (forced to `True` if `markers` is `False`).",
    ],
    log_x=[
        "boolean (default `False`)",
        "If `True`, the x-axis is log-scaled in cartesian coordinates.",
    ],
    log_y=[
        "boolean (default `False`)",
        "If `True`, the y-axis is log-scaled in cartesian coordinates.",
    ],
    log_z=[
        "boolean (default `False`)",
        "If `True`, the z-axis is log-scaled in cartesian coordinates.",
    ],
    log_r=[
        "boolean (default `False`)",
        "If `True`, the radial axis is log-scaled in polar coordinates.",
    ],
    range_x=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the x-axis in cartesian coordinates.",
    ],
    range_y=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the y-axis in cartesian coordinates.",
    ],
    range_z=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the z-axis in cartesian coordinates.",
    ],
    range_color=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the continuous color scale.",
    ],
    range_r=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the radial axis in polar coordinates.",
    ],
    range_theta=[
        "list of two numbers",
        "If provided, overrides auto-scaling on the angular axis in polar coordinates.",
    ],
    title=["str", "The figure title."],
    template=[
        "str or dict or plotly.graph_objects.layout.Template instance",
        "The figure template name (must be a key in plotly.io.templates) or definition.",
    ],
    width=["int (default `None`)", "The figure width in pixels."],
    height=["int (default `None`)", "The figure height in pixels."],
    labels=[
        "dict with str keys and str values (default `{}`)",
        "By default, column names are used in the figure for axis titles, legend entries and hovers.",
        "This parameter allows this to be overridden.",
        "The keys of this dict should correspond to column names, and the values should correspond to the desired label to be displayed.",
    ],
    category_orders=[
        "dict with str keys and list of str values (default `{}`)",
        "By default, in Python 3.6+, the order of categorical values in axes, legends and facets depends on the order in which these values are first encountered in `data_frame` (and no order is guaranteed by default in Python below 3.6).",
        "This parameter is used to force a specific ordering of values per column.",
        "The keys of this dict should correspond to column names, and the values should be lists of strings corresponding to the specific display order desired.",
    ],
    marginal=[
        "str",
        "One of `'rug'`, `'box'`, `'violin'`, or `'histogram'`.",
        "If set, a subplot is drawn alongside the main plot, visualizing the distribution.",
    ],
    marginal_x=[
        "str",
        "One of `'rug'`, `'box'`, `'violin'`, or `'histogram'`.",
        "If set, a horizontal subplot is drawn above the main plot, visualizing the x-distribution.",
    ],
    marginal_y=[
        "str",
        "One of `'rug'`, `'box'`, `'violin'`, or `'histogram'`.",
        "If set, a vertical subplot is drawn to the right of the main plot, visualizing the y-distribution.",
    ],
    trendline=[
        "str",
        "One of `'ols'`, `'lowess'`, `'rolling'`, `'expanding'` or `'ewm'`.",
        "If `'ols'`, an Ordinary Least Squares regression line will be drawn for each discrete-color/symbol group.",
        "If `'lowess`', a Locally Weighted Scatterplot Smoothing line will be drawn for each discrete-color/symbol group.",
        "If `'rolling`', a Rolling (e.g. rolling average, rolling median) line will be drawn for each discrete-color/symbol group.",
        "If `'expanding`', an Expanding (e.g. expanding average, expanding sum) line will be drawn for each discrete-color/symbol group.",
        "If `'ewm`', an Exponentially Weighted Moment (e.g. exponentially-weighted moving average) line will be drawn for each discrete-color/symbol group.",
        "See the docstrings for the functions in `plotly.express.trendline_functions` for more details on these functions and how",
        "to configure them with the `trendline_options` argument.",
    ],
    trendline_options=[
        "dict",
        "Options passed as the first argument to the function from `plotly.express.trendline_functions` ",
        "named in the `trendline` argument.",
    ],
    trendline_color_override=[
        "str",
        "Valid CSS color.",
        "If provided, and if `trendline` is set, all trendlines will be drawn in this color rather than in the same color as the traces from which they draw their inputs.",
    ],
    trendline_scope=[
        "str (one of `'trace'` or `'overall'`, default `'trace'`)",
        "If `'trace'`, then one trendline is drawn per trace (i.e. per color, symbol, facet, animation frame etc) and if `'overall'` then one trendline is computed for the entire dataset, and replicated across all facets.",
    ],
    render_mode=[
        "str",
        "One of `'auto'`, `'svg'` or `'webgl'`, default `'auto'`",
        "Controls the browser API used to draw marks.",
        "`'svg'` is appropriate for figures of less than 1000 data points, and will allow for fully-vectorized output.",
        "`'webgl'` is likely necessary for acceptable performance above 1000 points but rasterizes part of the output. ",
        "`'auto'` uses heuristics to choose the mode.",
    ],
    direction=[
        "str",
        "One of '`counterclockwise'` or `'clockwise'`. Default is `'clockwise'`",
        "Sets the direction in which increasing values of the angular axis are drawn.",
    ],
    start_angle=[
        "int (default `90`)",
        "Sets start angle for the angular axis, with 0 being due east and 90 being due north.",
    ],
    histfunc=[
        "str (default `'count'` if no arguments are provided, else `'sum'`)",
        "One of `'count'`, `'sum'`, `'avg'`, `'min'`, or `'max'`.",
        "Function used to aggregate values for summarization (note: can be normalized with `histnorm`).",
    ],
    histnorm=[
        "str (default `None`)",
        "One of `'percent'`, `'probability'`, `'density'`, or `'probability density'`",
        "If `None`, the output of `histfunc` is used as is.",
        "If `'probability'`, the output of `histfunc` for a given bin is divided by the sum of the output of `histfunc` for all bins.",
        "If `'percent'`, the output of `histfunc` for a given bin is divided by the sum of the output of `histfunc` for all bins and multiplied by 100.",
        "If `'density'`, the output of `histfunc` for a given bin is divided by the size of the bin.",
        "If `'probability density'`, the output of `histfunc` for a given bin is normalized such that it corresponds to the probability that a random event whose distribution is described by the output of `histfunc` will fall into that bin.",
    ],
    barnorm=[
        "str (default `None`)",
        "One of `'fraction'` or `'percent'`.",
        "If `'fraction'`, the value of each bar is divided by the sum of all values at that location coordinate.",
        "`'percent'` is the same but multiplied by 100 to show percentages.",
        "`None` will stack up all values at each location coordinate.",
    ],
    groupnorm=[
        "str (default `None`)",
        "One of `'fraction'` or `'percent'`.",
        "If `'fraction'`, the value of each point is divided by the sum of all values at that location coordinate.",
        "`'percent'` is the same but multiplied by 100 to show percentages.",
        "`None` will stack up all values at each location coordinate.",
    ],
    barmode=[
        "str (default `'relative'`)",
        "One of `'group'`, `'overlay'` or `'relative'`",
        "In `'relative'` mode, bars are stacked above zero for positive values and below zero for negative values.",
        "In `'overlay'` mode, bars are drawn on top of one another.",
        "In `'group'` mode, bars are placed beside each other.",
    ],
    boxmode=[
        "str (default `'group'`)",
        "One of `'group'` or `'overlay'`",
        "In `'overlay'` mode, boxes are on drawn top of one another.",
        "In `'group'` mode, boxes are placed beside each other.",
    ],
    violinmode=[
        "str (default `'group'`)",
        "One of `'group'` or `'overlay'`",
        "In `'overlay'` mode, violins are on drawn top of one another.",
        "In `'group'` mode, violins are placed beside each other.",
    ],
    stripmode=[
        "str (default `'group'`)",
        "One of `'group'` or `'overlay'`",
        "In `'overlay'` mode, strips are on drawn top of one another.",
        "In `'group'` mode, strips are placed beside each other.",
    ],
    zoom=["int (default `8`)", "Between 0 and 20.", "Sets map zoom level."],
    orientation=[
        "str, one of `'h'` for horizontal or `'v'` for vertical. ",
        "(default `'v'` if `x` and `y` are provided and both continous or both categorical, ",
        "otherwise `'v'`(`'h'`) if `x`(`y`) is categorical and `y`(`x`) is continuous, ",
        "otherwise `'v'`(`'h'`) if only `x`(`y`) is provided) ",
    ],
    line_close=[
        "boolean (default `False`)",
        "If `True`, an extra line segment is drawn between the first and last point.",
    ],
    line_shape=[
        "str (default `'linear'`)",
        "One of `'linear'`, `'spline'`, `'hv'`, `'vh'`, `'hvh'`, or `'vhv'`",
    ],
    fitbounds=["str (default `False`).", "One of `False`, `locations` or `geojson`."],
    basemap_visible=["bool", "Force the basemap visibility."],
    scope=[
        "str (default `'world'`).",
        "One of `'world'`, `'usa'`, `'europe'`, `'asia'`, `'africa'`, `'north america'`, or `'south america'`"
        "Default is `'world'` unless `projection` is set to `'albers usa'`, which forces `'usa'`.",
    ],
    projection=[
        "str ",
        "One of `'equirectangular'`, `'mercator'`, `'orthographic'`, `'natural earth'`, `'kavrayskiy7'`, `'miller'`, `'robinson'`, `'eckert4'`, `'azimuthal equal area'`, `'azimuthal equidistant'`, `'conic equal area'`, `'conic conformal'`, `'conic equidistant'`, `'gnomonic'`, `'stereographic'`, `'mollweide'`, `'hammer'`, `'transverse mercator'`, `'albers usa'`, `'winkel tripel'`, `'aitoff'`, or `'sinusoidal'`"
        "Default depends on `scope`.",
    ],
    center=[
        "dict",
        "Dict keys are `'lat'` and `'lon'`",
        "Sets the center point of the map.",
    ],
    map_style=[
        "str (default `'basic'`)",
        "Identifier of base map style.",
        "Allowed values are `'basic'`, `'carto-darkmatter'`, `'carto-darkmatter-nolabels'`, `'carto-positron'`, `'carto-positron-nolabels'`, `'carto-voyager'`, `'carto-voyager-nolabels'`, `'dark'`, `'light'`, `'open-street-map'`, `'outdoors'`, `'satellite'`, `'satellite-streets'`, `'streets'`, `'white-bg'`.",
    ],
    mapbox_style=[
        "str (default `'basic'`, needs Mapbox API token)",
        "Identifier of base map style, some of which require a Mapbox or Stadia Maps API token to be set using `plotly.express.set_mapbox_access_token()`.",
        "Allowed values which do not require a token are `'open-street-map'`, `'white-bg'`, `'carto-positron'`, `'carto-darkmatter'`.",
        "Allowed values which require a Mapbox API token are `'basic'`, `'streets'`, `'outdoors'`, `'light'`, `'dark'`, `'satellite'`, `'satellite-streets'`.",
        "Allowed values which require a Stadia Maps API token are `'stamen-terrain'`, `'stamen-toner'`, `'stamen-watercolor'`.",
    ],
    points=[
        "str or boolean (default `'outliers'`)",
        "One of `'outliers'`, `'suspectedoutliers'`, `'all'`, or `False`.",
        "If `'outliers'`, only the sample points lying outside the whiskers are shown.",
        "If `'suspectedoutliers'`, all outlier points are shown and those less than 4*Q1-3*Q3 or greater than 4*Q3-3*Q1 are highlighted with the marker's `'outliercolor'`.",
        "If `'outliers'`, only the sample points lying outside the whiskers are shown.",
        "If `'all'`, all sample points are shown.",
        "If `False`, no sample points are shown and the whiskers extend to the full range of the sample.",
    ],
    box=["boolean (default `False`)", "If `True`, boxes are drawn inside the violins."],
    notched=["boolean (default `False`)", "If `True`, boxes are drawn with notches."],
    geojson=[
        "GeoJSON-formatted dict",
        "Must contain a Polygon feature collection, with IDs, which are references from `locations`.",
    ],
    featureidkey=[
        "str (default: `'id'`)",
        "Path to field in GeoJSON feature object with which to match the values passed in to `locations`."
        "The most common alternative to the default is of the form `'properties.<key>`.",
    ],
    cumulative=[
        "boolean (default `False`)",
        "If `True`, histogram values are cumulative.",
    ],
    nbins=["int", "Positive integer.", "Sets the number of bins."],
    nbinsx=["int", "Positive integer.", "Sets the number of bins along the x axis."],
    nbinsy=["int", "Positive integer.", "Sets the number of bins along the y axis."],
    branchvalues=[
        "str",
        "'total' or 'remainder'",
        "Determines how the items in `values` are summed. When"
        "set to 'total', items in `values` are taken to be value"
        "of all its descendants. When set to 'remainder', items"
        "in `values` corresponding to the root and the branches"
        ":sectors are taken to be the extra part not part of the"
        "sum of the values at their leaves.",
    ],
    maxdepth=[
        "int",
        "Positive integer",
        "Sets the number of rendered sectors from any given `level`. Set `maxdepth` to -1 to render all the"
        "levels in the hierarchy.",
    ],
    ecdfnorm=[
        "string or `None` (default `'probability'`)",
        "One of `'probability'` or `'percent'`",
        "If `None`, values will be raw counts or sums.",
        "If `'probability', values will be probabilities normalized from 0 to 1.",
        "If `'percent', values will be percentages normalized from 0 to 100.",
    ],
    ecdfmode=[
        "string (default `'standard'`)",
        "One of `'standard'`, `'complementary'` or `'reversed'`",
        "If `'standard'`, the ECDF is plotted such that values represent data at or below the point.",
        "If `'complementary'`, the CCDF is plotted such that values represent data above the point.",
        "If `'reversed'`, a variant of the CCDF is plotted such that values represent data at or above the point.",
    ],
    text_auto=[
        "bool or string (default `False`)",
        "If `True` or a string, the x or y or z values will be displayed as text, depending on the orientation",
        "A string like `'.2f'` will be interpreted as a `texttemplate` numeric formatting directive.",
    ],
)


def make_docstring(fn, override_dict=None, append_dict=None):
    override_dict = {} if override_dict is None else override_dict
    append_dict = {} if append_dict is None else append_dict
    tw = TextWrapper(width=75, initial_indent="    ", subsequent_indent="    ")
    result = (fn.__doc__ or "") + "\nParameters\n----------\n"
    for param in getfullargspec(fn)[0]:
        if override_dict.get(param):
            param_doc = list(override_dict[param])
        else:
            param_doc = list(docs[param])
            if append_dict.get(param):
                param_doc += append_dict[param]
        param_desc_list = param_doc[1:]
        param_desc = (
            tw.fill(" ".join(param_desc_list or ""))
            if param in docs or param in override_dict
            else "(documentation missing from map)"
        )

        param_type = param_doc[0]
        result += "%s: %s\n%s\n" % (param, param_type, param_desc)
    result += "\nReturns\n-------\n"
    result += "    plotly.graph_objects.Figure"
    return result
