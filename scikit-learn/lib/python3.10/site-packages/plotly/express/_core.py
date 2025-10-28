import plotly.graph_objs as go
import plotly.io as pio
from collections import namedtuple, OrderedDict
from ._special_inputs import IdentityMap, Constant, Range
from .trendline_functions import ols, lowess, rolling, expanding, ewm

from _plotly_utils.basevalidators import ColorscaleValidator
from plotly.colors import qualitative, sequential
import math

from plotly._subplots import (
    make_subplots,
    _set_trace_grid_reference,
    _subplot_type_for_trace_type,
)

import narwhals.stable.v1 as nw

# The reason to use narwhals.stable.v1 is to have a stable and perfectly
# backwards-compatible API, hence the confidence to not pin the Narwhals version exactly,
# allowing for multiple major libraries to have Narwhals as a dependency without
# forbidding users to install them all together due to dependency conflicts.

NO_COLOR = "px_no_color_constant"


trendline_functions = dict(
    lowess=lowess, rolling=rolling, ewm=ewm, expanding=expanding, ols=ols
)

# Declare all supported attributes, across all plot types
direct_attrables = (
    ["base", "x", "y", "z", "a", "b", "c", "r", "theta", "size", "x_start", "x_end"]
    + ["hover_name", "text", "names", "values", "parents", "wide_cross"]
    + ["ids", "error_x", "error_x_minus", "error_y", "error_y_minus", "error_z"]
    + ["error_z_minus", "lat", "lon", "locations", "animation_group"]
)
array_attrables = ["dimensions", "custom_data", "hover_data", "path", "wide_variable"]
group_attrables = ["animation_frame", "facet_row", "facet_col", "line_group"]
renameable_group_attrables = [
    "color",  # renamed to marker.color or line.color in infer_config
    "symbol",  # renamed to marker.symbol in infer_config
    "line_dash",  # renamed to line.dash in infer_config
    "pattern_shape",  # renamed to marker.pattern.shape in infer_config
]
all_attrables = (
    direct_attrables + array_attrables + group_attrables + renameable_group_attrables
)

cartesians = [go.Scatter, go.Scattergl, go.Bar, go.Funnel, go.Box, go.Violin]
cartesians += [go.Histogram, go.Histogram2d, go.Histogram2dContour]


class PxDefaults(object):
    __slots__ = [
        "template",
        "width",
        "height",
        "color_discrete_sequence",
        "color_discrete_map",
        "color_continuous_scale",
        "symbol_sequence",
        "symbol_map",
        "line_dash_sequence",
        "line_dash_map",
        "pattern_shape_sequence",
        "pattern_shape_map",
        "size_max",
        "category_orders",
        "labels",
    ]

    def __init__(self):
        self.reset()

    def reset(self):
        self.template = None
        self.width = None
        self.height = None
        self.color_discrete_sequence = None
        self.color_discrete_map = {}
        self.color_continuous_scale = None
        self.symbol_sequence = None
        self.symbol_map = {}
        self.line_dash_sequence = None
        self.line_dash_map = {}
        self.pattern_shape_sequence = None
        self.pattern_shape_map = {}
        self.size_max = 20
        self.category_orders = {}
        self.labels = {}


defaults = PxDefaults()
del PxDefaults


MAPBOX_TOKEN = None


def set_mapbox_access_token(token):
    """
    Arguments:
        token: A Mapbox token to be used in `plotly.express.scatter_mapbox` and \
        `plotly.express.line_mapbox` figures. See \
        https://docs.mapbox.com/help/how-mapbox-works/access-tokens/ for more details
    """
    global MAPBOX_TOKEN
    MAPBOX_TOKEN = token


def get_trendline_results(fig):
    """
    Extracts fit statistics for trendlines (when applied to figures generated with
    the `trendline` argument set to `"ols"`).

    Arguments:
        fig: the output of a `plotly.express` charting call
    Returns:
        A `pandas.DataFrame` with a column "px_fit_results" containing the `statsmodels`
        results objects, along with columns identifying the subset of the data the
        trendline was fit on.
    """
    return fig._px_trendlines


Mapping = namedtuple(
    "Mapping",
    [
        "show_in_trace_name",
        "grouper",
        "val_map",
        "sequence",
        "updater",
        "variable",
        "facet",
    ],
)
TraceSpec = namedtuple("TraceSpec", ["constructor", "attrs", "trace_patch", "marginal"])


def get_label(args, column):
    try:
        return args["labels"][column]
    except Exception:
        return column


def invert_label(args, column):
    """Invert mapping.
    Find key corresponding to value column in dict args["labels"].
    Returns `column` if the value does not exist.
    """
    reversed_labels = {value: key for (key, value) in args["labels"].items()}
    try:
        return reversed_labels[column]
    except Exception:
        return column


def _is_continuous(df: nw.DataFrame, col_name: str) -> bool:
    if nw.dependencies.is_pandas_like_dataframe(df_native := df.to_native()):
        # fastpath for pandas: Narwhals' Series.dtype has a bit of overhead, as it
        # tries to distinguish between true "object" columns, and "string" columns
        # disguised as "object". But here, we deal with neither.
        return df_native[col_name].dtype.kind in "ifc"
    return df.get_column(col_name).dtype.is_numeric()


def _to_unix_epoch_seconds(s: nw.Series) -> nw.Series:
    dtype = s.dtype
    if dtype == nw.Date:
        return s.dt.timestamp("ms") / 1_000
    if dtype == nw.Datetime:
        if dtype.time_unit in ("s", "ms"):
            return s.dt.timestamp("ms") / 1_000
        elif dtype.time_unit == "us":
            return s.dt.timestamp("us") / 1_000_000
        elif dtype.time_unit == "ns":
            return s.dt.timestamp("ns") / 1_000_000_000
        else:
            msg = "Unexpected dtype, please report a bug"
            raise ValueError(msg)
    else:
        msg = f"Expected Date or Datetime, got {dtype}"
        raise TypeError(msg)


def _generate_temporary_column_name(n_bytes, columns) -> str:
    """Wraps of Narwhals generate_temporary_column_name to generate a token
    which is guaranteed to not be in columns, nor in [col + token for col in columns]
    """
    counter = 0
    while True:
        # This is guaranteed to not be in columns by Narwhals
        token = nw.generate_temporary_column_name(n_bytes, columns=columns)

        # Now check that it is not in the [col + token for col in columns] list
        if token not in {f"{c}{token}" for c in columns}:
            return token

        counter += 1
        if counter > 100:
            msg = (
                "Internal Error: Plotly was not able to generate a column name with "
                f"{n_bytes=} and not in {columns}.\n"
                "Please report this to "
                "https://github.com/plotly/plotly.py/issues/new and we will try to "
                "replicate and fix it."
            )
            raise AssertionError(msg)


def get_decorated_label(args, column, role):
    original_label = label = get_label(args, column)
    if "histfunc" in args and (
        (role == "z")
        or (role == "x" and "orientation" in args and args["orientation"] == "h")
        or (role == "y" and "orientation" in args and args["orientation"] == "v")
    ):
        histfunc = args["histfunc"] or "count"
        if histfunc != "count":
            label = "%s of %s" % (histfunc, label)
        else:
            label = "count"

        if "histnorm" in args and args["histnorm"] is not None:
            if label == "count":
                label = args["histnorm"]
            else:
                histnorm = args["histnorm"]
                if histfunc == "sum":
                    if histnorm == "probability":
                        label = "%s of %s" % ("fraction", label)
                    elif histnorm == "percent":
                        label = "%s of %s" % (histnorm, label)
                    else:
                        label = "%s weighted by %s" % (histnorm, original_label)
                elif histnorm == "probability":
                    label = "%s of sum of %s" % ("fraction", label)
                elif histnorm == "percent":
                    label = "%s of sum of %s" % ("percent", label)
                else:
                    label = "%s of %s" % (histnorm, label)

        if "barnorm" in args and args["barnorm"] is not None:
            label = "%s (normalized as %s)" % (label, args["barnorm"])

    return label


def make_mapping(args, variable):
    if variable == "line_group" or variable == "animation_frame":
        return Mapping(
            show_in_trace_name=False,
            grouper=args[variable],
            val_map={},
            sequence=[""],
            variable=variable,
            updater=(lambda trace, v: v),
            facet=None,
        )
    if variable == "facet_row" or variable == "facet_col":
        letter = "x" if variable == "facet_col" else "y"
        return Mapping(
            show_in_trace_name=False,
            variable=letter,
            grouper=args[variable],
            val_map={},
            sequence=[i for i in range(1, 1000)],
            updater=(lambda trace, v: v),
            facet="row" if variable == "facet_row" else "col",
        )
    (parent, variable, *other_variables) = variable.split(".")
    vprefix = variable
    arg_name = variable
    if variable == "color":
        vprefix = "color_discrete"
    if variable == "dash":
        arg_name = "line_dash"
        vprefix = "line_dash"
    if variable in ["pattern", "shape"]:
        arg_name = "pattern_shape"
        vprefix = "pattern_shape"
    if args[vprefix + "_map"] == "identity":
        val_map = IdentityMap()
    else:
        val_map = args[vprefix + "_map"].copy()
    return Mapping(
        show_in_trace_name=True,
        variable=variable,
        grouper=args[arg_name],
        val_map=val_map,
        sequence=args[vprefix + "_sequence"],
        updater=lambda trace, v: trace.update(
            {parent: {".".join([variable] + other_variables): v}}
        ),
        facet=None,
    )


def make_trace_kwargs(args, trace_spec, trace_data, mapping_labels, sizeref):
    """Populates a dict with arguments to update trace

    Parameters
    ----------
    args : dict
        args to be used for the trace
    trace_spec : NamedTuple
        which kind of trace to be used (has constructor, marginal etc.
        attributes)
    trace_data : pandas DataFrame
        data
    mapping_labels : dict
        to be used for hovertemplate
    sizeref : float
        marker sizeref

    Returns
    -------
    trace_patch : dict
        dict to be used to update trace
    fit_results : dict
        fit information to be used for trendlines
    """
    trace_data: nw.DataFrame
    df: nw.DataFrame = args["data_frame"]

    if "line_close" in args and args["line_close"]:
        trace_data = nw.concat([trace_data, trace_data.head(1)], how="vertical")

    trace_patch = trace_spec.trace_patch.copy() or {}
    fit_results = None
    hover_header = ""
    for attr_name in trace_spec.attrs:
        attr_value = args[attr_name]
        attr_label = get_decorated_label(args, attr_value, attr_name)
        if attr_name == "dimensions":
            dims = [
                (name, trace_data.get_column(name))
                for name in trace_data.columns
                if ((not attr_value) or (name in attr_value))
                and (trace_spec.constructor != go.Parcoords or _is_continuous(df, name))
                and (
                    trace_spec.constructor != go.Parcats
                    or (attr_value is not None and name in attr_value)
                    or nw.to_py_scalar(df.get_column(name).n_unique())
                    <= args["dimensions_max_cardinality"]
                )
            ]
            trace_patch["dimensions"] = [
                dict(label=get_label(args, name), values=column)
                for (name, column) in dims
            ]
            if trace_spec.constructor == go.Splom:
                for d in trace_patch["dimensions"]:
                    d["axis"] = dict(matches=True)
                mapping_labels["%{xaxis.title.text}"] = "%{x}"
                mapping_labels["%{yaxis.title.text}"] = "%{y}"

        elif attr_value is not None:
            if attr_name == "size":
                if "marker" not in trace_patch:
                    trace_patch["marker"] = dict()
                trace_patch["marker"]["size"] = trace_data.get_column(attr_value)
                trace_patch["marker"]["sizemode"] = "area"
                trace_patch["marker"]["sizeref"] = sizeref
                mapping_labels[attr_label] = "%{marker.size}"
            elif attr_name == "marginal_x":
                if trace_spec.constructor == go.Histogram:
                    mapping_labels["count"] = "%{y}"
            elif attr_name == "marginal_y":
                if trace_spec.constructor == go.Histogram:
                    mapping_labels["count"] = "%{x}"
            elif attr_name == "trendline":
                if (
                    args["x"]
                    and args["y"]
                    and len(
                        trace_data.select(nw.col(args["x"], args["y"])).drop_nulls()
                    )
                    > 1
                ):
                    # sorting is bad but trace_specs with "trendline" have no other attrs
                    sorted_trace_data = trace_data.sort(by=args["x"], nulls_last=True)
                    y = sorted_trace_data.get_column(args["y"])
                    x = sorted_trace_data.get_column(args["x"])

                    if x.dtype == nw.Datetime or x.dtype == nw.Date:
                        # convert to unix epoch seconds
                        x = _to_unix_epoch_seconds(x)
                    elif not x.dtype.is_numeric():
                        try:
                            x = x.cast(nw.Float64())
                        except ValueError:
                            raise ValueError(
                                "Could not convert value of 'x' ('%s') into a numeric type. "
                                "If 'x' contains stringified dates, please convert to a datetime column."
                                % args["x"]
                            )

                    if not y.dtype.is_numeric():
                        try:
                            y = y.cast(nw.Float64())
                        except ValueError:
                            raise ValueError(
                                "Could not convert value of 'y' into a numeric type."
                            )

                    # preserve original values of "x" in case they're dates
                    # otherwise numpy/pandas can mess with the timezones
                    # NB this means trendline functions must output one-to-one with the input series
                    # i.e. we can't do resampling, because then the X values might not line up!
                    non_missing = ~(x.is_null() | y.is_null())
                    trace_patch["x"] = sorted_trace_data.filter(non_missing).get_column(
                        args["x"]
                    )
                    if (
                        trace_patch["x"].dtype == nw.Datetime
                        and trace_patch["x"].dtype.time_zone is not None
                    ):
                        # Remove time zone so that local time is displayed
                        trace_patch["x"] = (
                            trace_patch["x"].dt.replace_time_zone(None).to_numpy()
                        )
                    else:
                        trace_patch["x"] = trace_patch["x"].to_numpy()

                    trendline_function = trendline_functions[attr_value]
                    y_out, hover_header, fit_results = trendline_function(
                        args["trendline_options"],
                        sorted_trace_data.get_column(args["x"]),  # narwhals series
                        x.to_numpy(),  # numpy array
                        y.to_numpy(),  # numpy array
                        args["x"],
                        args["y"],
                        non_missing.to_numpy(),  # numpy array
                    )
                    assert len(y_out) == len(trace_patch["x"]), (
                        "missing-data-handling failure in trendline code"
                    )
                    trace_patch["y"] = y_out
                    mapping_labels[get_label(args, args["x"])] = "%{x}"
                    mapping_labels[get_label(args, args["y"])] = "%{y} <b>(trend)</b>"
            elif attr_name.startswith("error"):
                error_xy = attr_name[:7]
                arr = "arrayminus" if attr_name.endswith("minus") else "array"
                if error_xy not in trace_patch:
                    trace_patch[error_xy] = {}
                trace_patch[error_xy][arr] = trace_data.get_column(attr_value)
            elif attr_name == "custom_data":
                if len(attr_value) > 0:
                    # here we store a data frame in customdata, and it's serialized
                    # as a list of row lists, which is what we want
                    trace_patch["customdata"] = trace_data.select(nw.col(attr_value))
            elif attr_name == "hover_name":
                if trace_spec.constructor not in [
                    go.Histogram,
                    go.Histogram2d,
                    go.Histogram2dContour,
                ]:
                    trace_patch["hovertext"] = trace_data.get_column(attr_value)
                    if hover_header == "":
                        hover_header = "<b>%{hovertext}</b><br><br>"
            elif attr_name == "hover_data":
                if trace_spec.constructor not in [
                    go.Histogram,
                    go.Histogram2d,
                    go.Histogram2dContour,
                ]:
                    hover_is_dict = isinstance(attr_value, dict)
                    customdata_cols = args.get("custom_data") or []
                    for col in attr_value:
                        if hover_is_dict and not attr_value[col]:
                            continue
                        if col in [
                            args.get("x"),
                            args.get("y"),
                            args.get("z"),
                            args.get("base"),
                        ]:
                            continue
                        try:
                            position = args["custom_data"].index(col)
                        except (ValueError, AttributeError, KeyError):
                            position = len(customdata_cols)
                            customdata_cols.append(col)
                        attr_label_col = get_decorated_label(args, col, None)
                        mapping_labels[attr_label_col] = "%%{customdata[%d]}" % (
                            position
                        )

                    if len(customdata_cols) > 0:
                        # here we store a data frame in customdata, and it's serialized
                        # as a list of row lists, which is what we want

                        # dict.fromkeys(customdata_cols) allows to deduplicate column
                        # names, yet maintaining the original order.
                        trace_patch["customdata"] = trace_data.select(
                            *[nw.col(c) for c in dict.fromkeys(customdata_cols)]
                        )
            elif attr_name == "color":
                if trace_spec.constructor in [
                    go.Choropleth,
                    go.Choroplethmap,
                    go.Choroplethmapbox,
                ]:
                    trace_patch["z"] = trace_data.get_column(attr_value)
                    trace_patch["coloraxis"] = "coloraxis1"
                    mapping_labels[attr_label] = "%{z}"
                elif trace_spec.constructor in [
                    go.Sunburst,
                    go.Treemap,
                    go.Icicle,
                    go.Pie,
                    go.Funnelarea,
                ]:
                    if "marker" not in trace_patch:
                        trace_patch["marker"] = dict()

                    if args.get("color_is_continuous"):
                        trace_patch["marker"]["colors"] = trace_data.get_column(
                            attr_value
                        )
                        trace_patch["marker"]["coloraxis"] = "coloraxis1"
                        mapping_labels[attr_label] = "%{color}"
                    else:
                        trace_patch["marker"]["colors"] = []
                        if args["color_discrete_map"] is not None:
                            mapping = args["color_discrete_map"].copy()
                        else:
                            mapping = {}
                        for cat in trace_data.get_column(attr_value).to_list():
                            # although trace_data.get_column(attr_value) is a Narwhals
                            # Series, which is an iterable, explicitly calling a to_list()
                            # makes sure that the elements we loop over are python objects
                            # in all cases, since depending on the backend this may not be
                            # the case (e.g. PyArrow)
                            if mapping.get(cat) is None:
                                mapping[cat] = args["color_discrete_sequence"][
                                    len(mapping) % len(args["color_discrete_sequence"])
                                ]
                            trace_patch["marker"]["colors"].append(mapping[cat])
                else:
                    colorable = "marker"
                    if trace_spec.constructor in [go.Parcats, go.Parcoords]:
                        colorable = "line"
                    if colorable not in trace_patch:
                        trace_patch[colorable] = dict()
                    trace_patch[colorable]["color"] = trace_data.get_column(attr_value)
                    trace_patch[colorable]["coloraxis"] = "coloraxis1"
                    mapping_labels[attr_label] = "%%{%s.color}" % colorable
            elif attr_name == "animation_group":
                trace_patch["ids"] = trace_data.get_column(attr_value)
            elif attr_name == "locations":
                trace_patch[attr_name] = trace_data.get_column(attr_value)
                mapping_labels[attr_label] = "%{location}"
            elif attr_name == "values":
                trace_patch[attr_name] = trace_data.get_column(attr_value)
                _label = "value" if attr_label == "values" else attr_label
                mapping_labels[_label] = "%{value}"
            elif attr_name == "parents":
                trace_patch[attr_name] = trace_data.get_column(attr_value)
                _label = "parent" if attr_label == "parents" else attr_label
                mapping_labels[_label] = "%{parent}"
            elif attr_name == "ids":
                trace_patch[attr_name] = trace_data.get_column(attr_value)
                _label = "id" if attr_label == "ids" else attr_label
                mapping_labels[_label] = "%{id}"
            elif attr_name == "names":
                if trace_spec.constructor in [
                    go.Sunburst,
                    go.Treemap,
                    go.Icicle,
                    go.Pie,
                    go.Funnelarea,
                ]:
                    trace_patch["labels"] = trace_data.get_column(attr_value)
                    _label = "label" if attr_label == "names" else attr_label
                    mapping_labels[_label] = "%{label}"
                else:
                    trace_patch[attr_name] = trace_data.get_column(attr_value)
            else:
                trace_patch[attr_name] = trace_data.get_column(attr_value)
                mapping_labels[attr_label] = "%%{%s}" % attr_name
        elif (trace_spec.constructor == go.Histogram and attr_name in ["x", "y"]) or (
            trace_spec.constructor in [go.Histogram2d, go.Histogram2dContour]
            and attr_name == "z"
        ):
            # ensure that stuff like "count" gets into the hoverlabel
            mapping_labels[attr_label] = "%%{%s}" % attr_name
    if trace_spec.constructor not in [go.Parcoords, go.Parcats]:
        # Modify mapping_labels according to hover_data keys
        # if hover_data is a dict
        mapping_labels_copy = OrderedDict(mapping_labels)
        if args["hover_data"] and isinstance(args["hover_data"], dict):
            for k, v in mapping_labels.items():
                # We need to invert the mapping here
                k_args = invert_label(args, k)
                if k_args in args["hover_data"]:
                    formatter = args["hover_data"][k_args][0]
                    if formatter:
                        if isinstance(formatter, str):
                            mapping_labels_copy[k] = v.replace("}", "%s}" % formatter)
                    else:
                        _ = mapping_labels_copy.pop(k)
        hover_lines = [k + "=" + v for k, v in mapping_labels_copy.items()]
        trace_patch["hovertemplate"] = hover_header + "<br>".join(hover_lines)
        trace_patch["hovertemplate"] += "<extra></extra>"
    return trace_patch, fit_results


def configure_axes(args, constructor, fig, orders):
    configurators = {
        go.Scatter3d: configure_3d_axes,
        go.Scatterternary: configure_ternary_axes,
        go.Scatterpolar: configure_polar_axes,
        go.Scatterpolargl: configure_polar_axes,
        go.Barpolar: configure_polar_axes,
        go.Scattermap: configure_map,
        go.Choroplethmap: configure_map,
        go.Densitymap: configure_map,
        go.Scattermapbox: configure_mapbox,
        go.Choroplethmapbox: configure_mapbox,
        go.Densitymapbox: configure_mapbox,
        go.Scattergeo: configure_geo,
        go.Choropleth: configure_geo,
    }
    for c in cartesians:
        configurators[c] = configure_cartesian_axes
    if constructor in configurators:
        configurators[constructor](args, fig, orders)


def set_cartesian_axis_opts(args, axis, letter, orders):
    log_key = "log_" + letter
    range_key = "range_" + letter
    if log_key in args and args[log_key]:
        axis["type"] = "log"
        if range_key in args and args[range_key]:
            axis["range"] = [math.log(r, 10) for r in args[range_key]]
    elif range_key in args and args[range_key]:
        axis["range"] = args[range_key]

    if args[letter] in orders:
        axis["categoryorder"] = "array"
        axis["categoryarray"] = (
            orders[args[letter]]
            if isinstance(axis, go.layout.XAxis)
            else list(reversed(orders[args[letter]]))  # top down for Y axis
        )


def configure_cartesian_marginal_axes(args, fig, orders):
    nrows = len(fig._grid_ref)
    ncols = len(fig._grid_ref[0])

    # Set y-axis titles and axis options in the left-most column
    for yaxis in fig.select_yaxes(col=1):
        set_cartesian_axis_opts(args, yaxis, "y", orders)

    # Set x-axis titles and axis options in the bottom-most row
    for xaxis in fig.select_xaxes(row=1):
        set_cartesian_axis_opts(args, xaxis, "x", orders)

    # Configure axis ticks on marginal subplots
    if args["marginal_x"]:
        fig.update_yaxes(
            showticklabels=False, showline=False, ticks="", range=None, row=nrows
        )
        if args["template"].layout.yaxis.showgrid is None:
            fig.update_yaxes(showgrid=args["marginal_x"] == "histogram", row=nrows)
        if args["template"].layout.xaxis.showgrid is None:
            fig.update_xaxes(showgrid=True, row=nrows)

    if args["marginal_y"]:
        fig.update_xaxes(
            showticklabels=False, showline=False, ticks="", range=None, col=ncols
        )
        if args["template"].layout.xaxis.showgrid is None:
            fig.update_xaxes(showgrid=args["marginal_y"] == "histogram", col=ncols)
        if args["template"].layout.yaxis.showgrid is None:
            fig.update_yaxes(showgrid=True, col=ncols)

    # Add axis titles to non-marginal subplots
    y_title = get_decorated_label(args, args["y"], "y")
    if args["marginal_x"]:
        fig.update_yaxes(title_text=y_title, row=1, col=1)
    else:
        for row in range(1, nrows + 1):
            fig.update_yaxes(title_text=y_title, row=row, col=1)

    x_title = get_decorated_label(args, args["x"], "x")
    if args["marginal_y"]:
        fig.update_xaxes(title_text=x_title, row=1, col=1)
    else:
        for col in range(1, ncols + 1):
            fig.update_xaxes(title_text=x_title, row=1, col=col)

    # Configure axis type across all x-axes
    if "log_x" in args and args["log_x"]:
        fig.update_xaxes(type="log")

    # Configure axis type across all y-axes
    if "log_y" in args and args["log_y"]:
        fig.update_yaxes(type="log")

    # Configure matching and axis type for marginal y-axes
    matches_y = "y" + str(ncols + 1)
    if args["marginal_x"]:
        for row in range(2, nrows + 1, 2):
            fig.update_yaxes(matches=matches_y, type=None, row=row)

    if args["marginal_y"]:
        for col in range(2, ncols + 1, 2):
            fig.update_xaxes(matches="x2", type=None, col=col)


def configure_cartesian_axes(args, fig, orders):
    if ("marginal_x" in args and args["marginal_x"]) or (
        "marginal_y" in args and args["marginal_y"]
    ):
        configure_cartesian_marginal_axes(args, fig, orders)
        return

    # Set y-axis titles and axis options in the left-most column
    y_title = get_decorated_label(args, args["y"], "y")
    for yaxis in fig.select_yaxes(col=1):
        yaxis.update(title_text=y_title)
        set_cartesian_axis_opts(args, yaxis, "y", orders)

    # Set x-axis titles and axis options in the bottom-most row
    x_title = get_decorated_label(args, args["x"], "x")
    for xaxis in fig.select_xaxes(row=1):
        if "is_timeline" not in args:
            xaxis.update(title_text=x_title)
        set_cartesian_axis_opts(args, xaxis, "x", orders)

    # Configure axis type across all x-axes
    if "log_x" in args and args["log_x"]:
        fig.update_xaxes(type="log")

    # Configure axis type across all y-axes
    if "log_y" in args and args["log_y"]:
        fig.update_yaxes(type="log")

    if "is_timeline" in args:
        fig.update_xaxes(type="date")

    if "ecdfmode" in args:
        if args["orientation"] == "v":
            fig.update_yaxes(rangemode="tozero")
        else:
            fig.update_xaxes(rangemode="tozero")


def configure_ternary_axes(args, fig, orders):
    fig.update_ternaries(
        aaxis=dict(title_text=get_label(args, args["a"])),
        baxis=dict(title_text=get_label(args, args["b"])),
        caxis=dict(title_text=get_label(args, args["c"])),
    )


def configure_polar_axes(args, fig, orders):
    patch = dict(
        angularaxis=dict(direction=args["direction"], rotation=args["start_angle"]),
        radialaxis=dict(),
    )

    for var, axis in [("r", "radialaxis"), ("theta", "angularaxis")]:
        if args[var] in orders:
            patch[axis]["categoryorder"] = "array"
            patch[axis]["categoryarray"] = orders[args[var]]

    radialaxis = patch["radialaxis"]
    if args["log_r"]:
        radialaxis["type"] = "log"
        if args["range_r"]:
            radialaxis["range"] = [math.log(x, 10) for x in args["range_r"]]
    else:
        if args["range_r"]:
            radialaxis["range"] = args["range_r"]

    if args["range_theta"]:
        patch["sector"] = args["range_theta"]
    fig.update_polars(patch)


def configure_3d_axes(args, fig, orders):
    patch = dict(
        xaxis=dict(title_text=get_label(args, args["x"])),
        yaxis=dict(title_text=get_label(args, args["y"])),
        zaxis=dict(title_text=get_label(args, args["z"])),
    )

    for letter in ["x", "y", "z"]:
        axis = patch[letter + "axis"]
        if args["log_" + letter]:
            axis["type"] = "log"
            if args["range_" + letter]:
                axis["range"] = [math.log(x, 10) for x in args["range_" + letter]]
        else:
            if args["range_" + letter]:
                axis["range"] = args["range_" + letter]
        if args[letter] in orders:
            axis["categoryorder"] = "array"
            axis["categoryarray"] = orders[args[letter]]
    fig.update_scenes(patch)


def configure_mapbox(args, fig, orders):
    center = args["center"]
    if not center and "lat" in args and "lon" in args:
        center = dict(
            lat=args["data_frame"][args["lat"]].mean(),
            lon=args["data_frame"][args["lon"]].mean(),
        )
    fig.update_mapboxes(
        accesstoken=MAPBOX_TOKEN,
        center=center,
        zoom=args["zoom"],
        style=args["mapbox_style"],
    )


def configure_map(args, fig, orders):
    center = args["center"]
    if not center and "lat" in args and "lon" in args:
        center = dict(
            lat=args["data_frame"][args["lat"]].mean(),
            lon=args["data_frame"][args["lon"]].mean(),
        )
    fig.update_maps(
        center=center,
        zoom=args["zoom"],
        style=args["map_style"],
    )


def configure_geo(args, fig, orders):
    fig.update_geos(
        center=args["center"],
        scope=args["scope"],
        fitbounds=args["fitbounds"],
        visible=args["basemap_visible"],
        projection=dict(type=args["projection"]),
    )


def configure_animation_controls(args, constructor, fig):
    def frame_args(duration):
        return {
            "frame": {"duration": duration, "redraw": constructor != go.Scatter},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }

    if "animation_frame" in args and args["animation_frame"] and len(fig.frames) > 1:
        fig.layout.updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, frame_args(500)],
                        "label": "&#9654;",
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ]
        fig.layout.sliders = [
            {
                "active": 0,
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {
                    "prefix": get_label(args, args["animation_frame"]) + "="
                },
                "pad": {"b": 10, "t": 60},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [[f.name], frame_args(0)],
                        "label": f.name,
                        "method": "animate",
                    }
                    for f in fig.frames
                ],
            }
        ]


def make_trace_spec(args, constructor, attrs, trace_patch):
    if constructor in [go.Scatter, go.Scatterpolar]:
        if "render_mode" in args and (
            args["render_mode"] == "webgl"
            or (
                args["render_mode"] == "auto"
                and len(args["data_frame"]) > 1000
                and args.get("line_shape") != "spline"
                and args["animation_frame"] is None
            )
        ):
            if constructor == go.Scatter:
                constructor = go.Scattergl
                if "orientation" in trace_patch:
                    del trace_patch["orientation"]
            else:
                constructor = go.Scatterpolargl
    # Create base trace specification
    result = [TraceSpec(constructor, attrs, trace_patch, None)]

    # Add marginal trace specifications
    for letter in ["x", "y"]:
        if "marginal_" + letter in args and args["marginal_" + letter]:
            trace_spec = None
            axis_map = dict(
                xaxis="x1" if letter == "x" else "x2",
                yaxis="y1" if letter == "y" else "y2",
            )
            if args["marginal_" + letter] == "histogram":
                trace_spec = TraceSpec(
                    constructor=go.Histogram,
                    attrs=[letter, "marginal_" + letter],
                    trace_patch=dict(opacity=0.5, bingroup=letter, **axis_map),
                    marginal=letter,
                )
            elif args["marginal_" + letter] == "violin":
                trace_spec = TraceSpec(
                    constructor=go.Violin,
                    attrs=[letter, "hover_name", "hover_data"],
                    trace_patch=dict(scalegroup=letter),
                    marginal=letter,
                )
            elif args["marginal_" + letter] == "box":
                trace_spec = TraceSpec(
                    constructor=go.Box,
                    attrs=[letter, "hover_name", "hover_data"],
                    trace_patch=dict(notched=True),
                    marginal=letter,
                )
            elif args["marginal_" + letter] == "rug":
                symbols = {"x": "line-ns-open", "y": "line-ew-open"}
                trace_spec = TraceSpec(
                    constructor=go.Box,
                    attrs=[letter, "hover_name", "hover_data"],
                    trace_patch=dict(
                        fillcolor="rgba(255,255,255,0)",
                        line={"color": "rgba(255,255,255,0)"},
                        boxpoints="all",
                        jitter=0,
                        hoveron="points",
                        marker={"symbol": symbols[letter]},
                    ),
                    marginal=letter,
                )
            if "color" in attrs or "color" not in args:
                if "marker" not in trace_spec.trace_patch:
                    trace_spec.trace_patch["marker"] = dict()
                first_default_color = args["color_continuous_scale"][0]
                trace_spec.trace_patch["marker"]["color"] = first_default_color
            result.append(trace_spec)

    # Add trendline trace specifications
    if args.get("trendline") and args.get("trendline_scope", "trace") == "trace":
        result.append(make_trendline_spec(args, constructor))
    return result


def make_trendline_spec(args, constructor):
    trace_spec = TraceSpec(
        constructor=(
            go.Scattergl
            if constructor == go.Scattergl  # could be contour
            else go.Scatter
        ),
        attrs=["trendline"],
        trace_patch=dict(mode="lines"),
        marginal=None,
    )
    if args["trendline_color_override"]:
        trace_spec.trace_patch["line"] = dict(color=args["trendline_color_override"])
    return trace_spec


def one_group(x):
    return ""


def apply_default_cascade(args):
    # first we apply px.defaults to unspecified args

    for param in defaults.__slots__:
        if param in args and args[param] is None:
            args[param] = getattr(defaults, param)

    # load the default template if set, otherwise "plotly"
    if args["template"] is None:
        if pio.templates.default is not None:
            args["template"] = pio.templates.default
        else:
            args["template"] = "plotly"

    try:
        # retrieve the actual template if we were given a name
        args["template"] = pio.templates[args["template"]]
    except Exception:
        # otherwise try to build a real template
        args["template"] = go.layout.Template(args["template"])

    # if colors not set explicitly or in px.defaults, defer to a template
    # if the template doesn't have one, we set some final fallback defaults
    if "color_continuous_scale" in args:
        if (
            args["color_continuous_scale"] is None
            and args["template"].layout.colorscale.sequential
        ):
            args["color_continuous_scale"] = [
                x[1] for x in args["template"].layout.colorscale.sequential
            ]
        if args["color_continuous_scale"] is None:
            args["color_continuous_scale"] = sequential.Viridis

    if "color_discrete_sequence" in args:
        if args["color_discrete_sequence"] is None and args["template"].layout.colorway:
            args["color_discrete_sequence"] = args["template"].layout.colorway
        if args["color_discrete_sequence"] is None:
            args["color_discrete_sequence"] = qualitative.D3

    # if symbol_sequence/line_dash_sequence not set explicitly or in px.defaults,
    # see if we can defer to template. If not, set reasonable defaults
    if "symbol_sequence" in args:
        if args["symbol_sequence"] is None and args["template"].data.scatter:
            args["symbol_sequence"] = [
                scatter.marker.symbol for scatter in args["template"].data.scatter
            ]
        if not args["symbol_sequence"] or not any(args["symbol_sequence"]):
            args["symbol_sequence"] = ["circle", "diamond", "square", "x", "cross"]

    if "line_dash_sequence" in args:
        if args["line_dash_sequence"] is None and args["template"].data.scatter:
            args["line_dash_sequence"] = [
                scatter.line.dash for scatter in args["template"].data.scatter
            ]
        if not args["line_dash_sequence"] or not any(args["line_dash_sequence"]):
            args["line_dash_sequence"] = [
                "solid",
                "dot",
                "dash",
                "longdash",
                "dashdot",
                "longdashdot",
            ]

    if "pattern_shape_sequence" in args:
        if args["pattern_shape_sequence"] is None and args["template"].data.bar:
            args["pattern_shape_sequence"] = [
                bar.marker.pattern.shape for bar in args["template"].data.bar
            ]
        if not args["pattern_shape_sequence"] or not any(
            args["pattern_shape_sequence"]
        ):
            args["pattern_shape_sequence"] = ["", "/", "\\", "x", "+", "."]


def _check_name_not_reserved(field_name, reserved_names):
    if field_name not in reserved_names:
        return field_name
    else:
        raise NameError(
            "A name conflict was encountered for argument '%s'. "
            "A column or index with name '%s' is ambiguous." % (field_name, field_name)
        )


def _get_reserved_col_names(args):
    """
    This function builds a list of columns of the data_frame argument used
    as arguments, either as str/int arguments or given as columns
    (pandas series type).
    """
    df: nw.DataFrame = args["data_frame"]
    reserved_names = set()
    for field in args:
        if field not in all_attrables:
            continue
        names = args[field] if field in array_attrables else [args[field]]
        if names is None:
            continue
        for arg in names:
            if arg is None:
                continue
            elif isinstance(arg, str):  # no need to add ints since kw arg are not ints
                reserved_names.add(arg)
            elif nw.dependencies.is_into_series(arg):
                arg_series = nw.from_native(arg, series_only=True)
                arg_name = arg_series.name
                if arg_name and arg_name in df.columns:
                    in_df = (arg_series == df.get_column(arg_name)).all()
                    if in_df:
                        reserved_names.add(arg_name)
            elif arg is nw.maybe_get_index(df) and arg.name is not None:
                reserved_names.add(arg.name)

    return reserved_names


def _is_col_list(columns, arg, is_pd_like, native_namespace):
    """Returns True if arg looks like it's a list of columns or references to columns
    in df_input, and False otherwise (in which case it's assumed to be a single column
    or reference to a column).
    """
    if arg is None or isinstance(arg, str) or isinstance(arg, int):
        return False
    if is_pd_like and isinstance(arg, native_namespace.MultiIndex):
        return False  # just to keep existing behaviour for now
    try:
        iter(arg)
    except TypeError:
        return False  # not iterable
    for c in arg:
        if isinstance(c, str) or isinstance(c, int):
            if columns is None or c not in columns:
                return False
        else:
            try:
                iter(c)
            except TypeError:
                return False  # not iterable
    return True


def _isinstance_listlike(x):
    """Returns True if x is an iterable which can be transformed into a pandas Series,
    False for the other types of possible values of a `hover_data` dict.
    A tuple of length 2 is a special case corresponding to a (format, data) tuple.
    """
    if (
        isinstance(x, str)
        or (isinstance(x, tuple) and len(x) == 2)
        or isinstance(x, bool)
        or x is None
    ):
        return False
    else:
        return True


def _escape_col_name(columns, col_name, extra):
    if columns is None:
        return col_name
    while col_name in columns or col_name in extra:
        col_name = "_" + col_name
    return col_name


def to_named_series(x, name=None, native_namespace=None):
    """Assuming x is list-like or even an existing Series, returns a new Series named `name`."""
    # With `pass_through=True`, the original object will be returned if unable to convert
    # to a Narwhals Series.
    x = nw.from_native(x, series_only=True, pass_through=True)
    if isinstance(x, nw.Series):
        return x.rename(name)
    elif native_namespace is not None:
        return nw.new_series(name=name, values=x, native_namespace=native_namespace)
    else:
        try:
            import pandas as pd

            return nw.new_series(name=name, values=x, native_namespace=pd)
        except ImportError:
            msg = "Pandas installation is required if no dataframe is provided."
            raise NotImplementedError(msg)


def process_args_into_dataframe(
    args, wide_mode, var_name, value_name, is_pd_like, native_namespace
):
    """
    After this function runs, the `all_attrables` keys of `args` all contain only
    references to columns of `df_output`. This function handles the extraction of data
    from `args["attrable"]` and column-name-generation as appropriate, and adds the
    data to `df_output` and then replaces `args["attrable"]` with the appropriate
    reference.
    """

    df_input: nw.DataFrame | None = args["data_frame"]
    df_provided = df_input is not None

    # we use a dict instead of a dataframe directly so that it doesn't cause
    # PerformanceWarning by pandas by repeatedly setting the columns.
    # a dict is used instead of a list as the columns needs to be overwritten.
    df_output = {}
    constants = {}
    ranges = []
    wide_id_vars = set()
    reserved_names = _get_reserved_col_names(args) if df_provided else set()

    # Case of functions with a "dimensions" kw: scatter_matrix, parcats, parcoords
    if "dimensions" in args and args["dimensions"] is None:
        if not df_provided:
            raise ValueError(
                "No data were provided. Please provide data either with the `data_frame` or with the `dimensions` argument."
            )
        else:
            df_output = {col: df_input.get_column(col) for col in df_input.columns}

    # hover_data is a dict
    hover_data_is_dict = (
        "hover_data" in args
        and args["hover_data"]
        and isinstance(args["hover_data"], dict)
    )
    # If dict, convert all values of hover_data to tuples to simplify processing
    if hover_data_is_dict:
        for k in args["hover_data"]:
            if _isinstance_listlike(args["hover_data"][k]):
                args["hover_data"][k] = (True, args["hover_data"][k])
            if not isinstance(args["hover_data"][k], tuple):
                args["hover_data"][k] = (args["hover_data"][k], None)
            if df_provided and args["hover_data"][k][1] is not None and k in df_input:
                raise ValueError(
                    "Ambiguous input: values for '%s' appear both in hover_data and data_frame"
                    % k
                )
    # Loop over possible arguments
    for field_name in all_attrables:
        # Massaging variables
        argument_list = (
            [args.get(field_name)]
            if field_name not in array_attrables
            else args.get(field_name)
        )

        # argument not specified, continue
        # The original also tested `or argument_list is [None]` but
        # that clause is always False, so it has been removed.  The
        # alternative fix would have been to test that `argument_list`
        # is of length 1 and its sole element is `None`, but that
        # feels pedantic. All tests pass with the change below; let's
        # see if the world decides we were wrong.
        if argument_list is None:
            continue

        # Argument name: field_name if the argument is not a list
        # Else we give names like ["hover_data_0, hover_data_1"] etc.
        field_list = (
            [field_name]
            if field_name not in array_attrables
            else [field_name + "_" + str(i) for i in range(len(argument_list))]
        )
        # argument_list and field_list ready, iterate over them
        # Core of the loop starts here
        for i, (argument, field) in enumerate(zip(argument_list, field_list)):
            length = len(df_output[next(iter(df_output))]) if len(df_output) else 0
            if argument is None:
                continue
            col_name = None
            # Case of multiindex
            if is_pd_like and isinstance(argument, native_namespace.MultiIndex):
                raise TypeError(
                    f"Argument '{field}' is a {native_namespace.__name__} MultiIndex. "
                    f"{native_namespace.__name__} MultiIndex is not supported by plotly "
                    "express at the moment."
                )
            # ----------------- argument is a special value ----------------------
            if isinstance(argument, (Constant, Range)):
                col_name = _check_name_not_reserved(
                    str(argument.label) if argument.label is not None else field,
                    reserved_names,
                )
                if isinstance(argument, Constant):
                    constants[col_name] = argument.value
                else:
                    ranges.append(col_name)
            # ----------------- argument is likely a col name ----------------------
            elif isinstance(argument, str) or not hasattr(argument, "__len__"):
                if (
                    field_name == "hover_data"
                    and hover_data_is_dict
                    and args["hover_data"][str(argument)][1] is not None
                ):
                    # hover_data has onboard data
                    # previously-checked to have no name-conflict with data_frame
                    col_name = str(argument)
                    real_argument = args["hover_data"][col_name][1]

                    if length and (real_length := len(real_argument)) != length:
                        raise ValueError(
                            "All arguments should have the same length. "
                            "The length of hover_data key `%s` is %d, whereas the "
                            "length of previously-processed arguments %s is %d"
                            % (
                                argument,
                                real_length,
                                str(list(df_output.keys())),
                                length,
                            )
                        )
                    df_output[col_name] = to_named_series(
                        real_argument, col_name, native_namespace
                    )
                elif not df_provided:
                    raise ValueError(
                        "String or int arguments are only possible when a "
                        "DataFrame or an array is provided in the `data_frame` "
                        "argument. No DataFrame was provided, but argument "
                        "'%s' is of type str or int." % field
                    )
                # Check validity of column name
                elif argument not in df_input.columns:
                    if wide_mode and argument in (value_name, var_name):
                        continue
                    else:
                        err_msg = (
                            "Value of '%s' is not the name of a column in 'data_frame'. "
                            "Expected one of %s but received: %s"
                            % (field, str(list(df_input.columns)), argument)
                        )
                        if argument == "index":
                            err_msg += "\n To use the index, pass it in directly as `df.index`."
                        raise ValueError(err_msg)
                elif length and (actual_len := len(df_input)) != length:
                    raise ValueError(
                        "All arguments should have the same length. "
                        "The length of column argument `df[%s]` is %d, whereas the "
                        "length of previously-processed arguments %s is %d"
                        % (
                            field,
                            actual_len,
                            str(list(df_output.keys())),
                            length,
                        )
                    )
                else:
                    col_name = str(argument)
                    df_output[col_name] = to_named_series(
                        df_input.get_column(argument), col_name
                    )
            # ----------------- argument is likely a column / array / list.... -------
            else:
                if df_provided and hasattr(argument, "name"):
                    if is_pd_like and argument is nw.maybe_get_index(df_input):
                        if argument.name is None or argument.name in df_input.columns:
                            col_name = "index"
                        else:
                            col_name = argument.name
                        col_name = _escape_col_name(
                            df_input.columns, col_name, [var_name, value_name]
                        )
                    else:
                        if (
                            argument.name is not None
                            and argument.name in df_input.columns
                            and (
                                to_named_series(
                                    argument, argument.name, native_namespace
                                )
                                == df_input.get_column(argument.name)
                            ).all()
                        ):
                            col_name = argument.name
                if col_name is None:  # numpy array, list...
                    col_name = _check_name_not_reserved(field, reserved_names)

                if length and (len_arg := len(argument)) != length:
                    raise ValueError(
                        "All arguments should have the same length. "
                        "The length of argument `%s` is %d, whereas the "
                        "length of previously-processed arguments %s is %d"
                        % (field, len_arg, str(list(df_output.keys())), length)
                    )

                df_output[str(col_name)] = to_named_series(
                    x=argument,
                    name=str(col_name),
                    native_namespace=native_namespace,
                )

            # Finally, update argument with column name now that column exists
            assert col_name is not None, (
                "Data-frame processing failure, likely due to a internal bug. "
                "Please report this to "
                "https://github.com/plotly/plotly.py/issues/new and we will try to "
                "replicate and fix it."
            )
            if field_name not in array_attrables:
                args[field_name] = str(col_name)
            elif isinstance(args[field_name], dict):
                pass
            else:
                args[field_name][i] = str(col_name)
            if field_name != "wide_variable":
                wide_id_vars.add(str(col_name))

    length = len(df_output[next(iter(df_output))]) if len(df_output) else 0

    if native_namespace is None:
        try:
            import pandas as pd

            native_namespace = pd
        except ImportError:
            msg = "Pandas installation is required if no dataframe is provided."
            raise NotImplementedError(msg)

    if ranges:
        import numpy as np

        range_series = nw.new_series(
            name="__placeholder__",
            values=np.arange(length),
            native_namespace=native_namespace,
        )
        df_output.update(
            {col_name: range_series.alias(col_name) for col_name in ranges}
        )

    df_output.update(
        {
            # constant is single value. repeat by len to avoid creating NaN on concatenating
            col_name: nw.new_series(
                name=col_name,
                values=[constants[col_name]] * length,
                native_namespace=native_namespace,
            )
            for col_name in constants
        }
    )

    if df_output:
        df_output = nw.from_dict(df_output)
    else:
        try:
            import pandas as pd
        except ImportError:
            msg = "Pandas installation is required."
            raise NotImplementedError(msg)
        df_output = nw.from_native(pd.DataFrame({}), eager_only=True)
    return df_output, wide_id_vars


def build_dataframe(args, constructor):
    """
    Constructs a dataframe and modifies `args` in-place.

    The argument values in `args` can be either strings corresponding to
    existing columns of a dataframe, or data arrays (lists, numpy arrays,
    pandas columns, series).

    Parameters
    ----------
    args : OrderedDict
        arguments passed to the px function and subsequently modified
    constructor : graph_object trace class
        the trace type selected for this figure
    """

    # make copies of all the fields via dict() and list()
    for field in args:
        if field in array_attrables and args[field] is not None:
            if isinstance(args[field], dict):
                args[field] = dict(args[field])
            elif field in ["custom_data", "hover_data"] and isinstance(
                args[field], str
            ):
                args[field] = [args[field]]
            else:
                args[field] = list(args[field])

    # Cast data_frame argument to DataFrame (it could be a numpy array, dict etc.)
    df_provided = args["data_frame"] is not None

    # Flag that indicates if the resulting data_frame after parsing is pandas-like
    # (in terms of resulting Narwhals DataFrame).
    # True if pandas, modin.pandas or cudf DataFrame/Series instance, or converted from
    # PySpark to pandas.
    is_pd_like = False

    # Flag that indicates if data_frame needs to be converted to PyArrow.
    # True if Ibis, DuckDB, Vaex, or implements __dataframe__
    needs_interchanging = False

    # If data_frame is provided, we parse it into a narwhals DataFrame, while accounting
    # for compatibility with pandas specific paths (e.g. Index/MultiIndex case).
    if df_provided:
        # data_frame is pandas-like DataFrame (pandas, modin.pandas, cudf)
        if nw.dependencies.is_pandas_like_dataframe(args["data_frame"]):
            columns = args["data_frame"].columns  # This can be multi index
            args["data_frame"] = nw.from_native(args["data_frame"], eager_only=True)
            is_pd_like = True

        # data_frame is pandas-like Series (pandas, modin.pandas, cudf)
        elif nw.dependencies.is_pandas_like_series(args["data_frame"]):
            args["data_frame"] = nw.from_native(
                args["data_frame"], series_only=True
            ).to_frame()
            columns = args["data_frame"].columns
            is_pd_like = True

        # data_frame is any other DataFrame object natively supported via Narwhals.
        # With `pass_through=True`, the original object will be returned if unable to convert
        # to a Narwhals DataFrame, making this condition False.
        elif isinstance(
            data_frame := nw.from_native(
                args["data_frame"], eager_or_interchange_only=True, pass_through=True
            ),
            nw.DataFrame,
        ):
            args["data_frame"] = data_frame
            needs_interchanging = nw.get_level(data_frame) == "interchange"
            columns = args["data_frame"].columns

        # data_frame is any other Series object natively supported via Narwhals.
        # With `pass_through=True`, the original object will be returned if unable to convert
        # to a Narwhals Series, making this condition False.
        elif isinstance(
            series := nw.from_native(
                args["data_frame"], series_only=True, pass_through=True
            ),
            nw.Series,
        ):
            args["data_frame"] = series.to_frame()
            columns = args["data_frame"].columns

        # data_frame is PySpark: it does not support interchange protocol and it is not
        # integrated in Narwhals. We use its native method to convert it to pandas.
        elif hasattr(args["data_frame"], "toPandas"):
            args["data_frame"] = nw.from_native(
                args["data_frame"].toPandas(), eager_only=True
            )
            columns = args["data_frame"].columns
            is_pd_like = True

        # data_frame is some other object type (e.g. dict, list, ...)
        # We try to import pandas, and then try to instantiate a pandas dataframe from
        # this such object
        else:
            try:
                import pandas as pd

                try:
                    args["data_frame"] = nw.from_native(
                        pd.DataFrame(args["data_frame"])
                    )
                    columns = args["data_frame"].columns
                    is_pd_like = True
                except Exception:
                    msg = (
                        f"Unable to convert data_frame of type {type(args['data_frame'])} "
                        "to pandas DataFrame. Please provide a supported dataframe type "
                        "or a type that can be passed to pd.DataFrame."
                    )

                    raise NotImplementedError(msg)
            except ImportError:
                msg = (
                    f"Attempting to convert data_frame of type {type(args['data_frame'])} "
                    "to pandas DataFrame, but Pandas is not installed. "
                    "Convert it to supported dataframe type or install pandas."
                )
                raise NotImplementedError(msg)

    # data_frame is not provided
    else:
        columns = None

    df_input: nw.DataFrame | None = args["data_frame"]
    index = (
        nw.maybe_get_index(df_input)
        if df_provided and not needs_interchanging
        else None
    )
    native_namespace = (
        nw.get_native_namespace(df_input)
        if df_provided and not needs_interchanging
        else None
    )

    # now we handle special cases like wide-mode or x-xor-y specification
    # by rearranging args to tee things up for process_args_into_dataframe to work
    no_x = args.get("x") is None
    no_y = args.get("y") is None
    wide_x = (
        False
        if no_x
        else _is_col_list(columns, args["x"], is_pd_like, native_namespace)
    )
    wide_y = (
        False
        if no_y
        else _is_col_list(columns, args["y"], is_pd_like, native_namespace)
    )

    wide_mode = False
    var_name = None  # will likely be "variable" in wide_mode
    wide_cross_name = None  # will likely be "index" in wide_mode
    value_name = None  # will likely be "value" in wide_mode
    hist2d_types = [go.Histogram2d, go.Histogram2dContour]
    hist1d_orientation = constructor == go.Histogram or "ecdfmode" in args
    if constructor in cartesians:
        if wide_x and wide_y:
            raise ValueError(
                "Cannot accept list of column references or list of columns for both `x` and `y`."
            )
        if df_provided and no_x and no_y:
            wide_mode = True
            if is_pd_like and isinstance(columns, native_namespace.MultiIndex):
                raise TypeError(
                    f"Data frame columns is a {native_namespace.__name__} MultiIndex. "
                    f"{native_namespace.__name__} MultiIndex is not supported by plotly "
                    "express at the moment."
                )
            args["wide_variable"] = list(columns)
            if is_pd_like and isinstance(columns, native_namespace.Index):
                var_name = columns.name
            else:
                var_name = None
            if var_name in [None, "value", "index"] or var_name in columns:
                var_name = "variable"
            if constructor == go.Funnel:
                wide_orientation = args.get("orientation") or "h"
            else:
                wide_orientation = args.get("orientation") or "v"
            args["orientation"] = wide_orientation
            args["wide_cross"] = None
        elif wide_x != wide_y:
            wide_mode = True
            args["wide_variable"] = args["y"] if wide_y else args["x"]
            if df_provided and is_pd_like and args["wide_variable"] is columns:
                var_name = columns.name
            if is_pd_like and isinstance(args["wide_variable"], native_namespace.Index):
                args["wide_variable"] = list(args["wide_variable"])
            if var_name in [None, "value", "index"] or (
                df_provided and var_name in columns
            ):
                var_name = "variable"
            if hist1d_orientation:
                wide_orientation = "v" if wide_x else "h"
            else:
                wide_orientation = "v" if wide_y else "h"
            args["y" if wide_y else "x"] = None
            args["wide_cross"] = None
            if not no_x and not no_y:
                wide_cross_name = "__x__" if wide_y else "__y__"

    if wide_mode:
        value_name = _escape_col_name(columns, "value", [])
        var_name = _escape_col_name(columns, var_name, [])

    # If the data_frame has interchange-only support levelin Narwhals, then we need to
    # convert it to a full support level backend.
    # Hence we convert requires Interchange to PyArrow.
    if needs_interchanging:
        if wide_mode:
            args["data_frame"] = nw.from_native(
                args["data_frame"].to_arrow(), eager_only=True
            )
        else:
            # Save precious resources by only interchanging columns that are
            # actually going to be plotted. This is tricky to do in the general case,
            # because Plotly allows calls like `px.line(df, x='x', y=['y1', df['y1']])`,
            # but interchange-only objects (e.g. DuckDB) don't typically have a concept
            # of self-standing Series. It's more important to perform project pushdown
            # here seeing as we're materialising to an (eager) PyArrow table.
            necessary_columns = {
                i for i in args.values() if isinstance(i, str) and i in columns
            }
            for field in args:
                if args[field] is not None and field in array_attrables:
                    necessary_columns.update(i for i in args[field] if i in columns)
            columns = list(necessary_columns)
            args["data_frame"] = nw.from_native(
                args["data_frame"].select(columns).to_arrow(), eager_only=True
            )
        import pyarrow as pa

        native_namespace = pa
    missing_bar_dim = None
    if (
        constructor in [go.Scatter, go.Bar, go.Funnel] + hist2d_types
        and not hist1d_orientation
    ):
        if not wide_mode and (no_x != no_y):
            for ax in ["x", "y"]:
                if args.get(ax) is None:
                    args[ax] = (
                        index
                        if index is not None
                        else Range(
                            label=_escape_col_name(columns, ax, [var_name, value_name])
                        )
                    )
                    if constructor == go.Bar:
                        missing_bar_dim = ax
                    else:
                        if args["orientation"] is None:
                            args["orientation"] = "v" if ax == "x" else "h"
        if wide_mode and wide_cross_name is None:
            if no_x != no_y and args["orientation"] is None:
                args["orientation"] = "v" if no_x else "h"
            if df_provided and is_pd_like and index is not None:
                if isinstance(index, native_namespace.MultiIndex):
                    raise TypeError(
                        f"Data frame index is a {native_namespace.__name__} MultiIndex. "
                        f"{native_namespace.__name__} MultiIndex is not supported by "
                        "plotly express at the moment."
                    )
                args["wide_cross"] = index
            else:
                args["wide_cross"] = Range(
                    label=_escape_col_name(columns, "index", [var_name, value_name])
                )

    no_color = False
    if isinstance(args.get("color"), str) and args["color"] == NO_COLOR:
        no_color = True
        args["color"] = None
    # now that things have been prepped, we do the systematic rewriting of `args`

    df_output, wide_id_vars = process_args_into_dataframe(
        args,
        wide_mode,
        var_name,
        value_name,
        is_pd_like,
        native_namespace,
    )
    df_output: nw.DataFrame
    # now that `df_output` exists and `args` contains only references, we complete
    # the special-case and wide-mode handling by further rewriting args and/or mutating
    # df_output

    count_name = _escape_col_name(df_output.columns, "count", [var_name, value_name])
    if not wide_mode and missing_bar_dim and constructor == go.Bar:
        # now that we've populated df_output, we check to see if the non-missing
        # dimension is categorical: if so, then setting the missing dimension to a
        # constant 1 is a less-insane thing to do than setting it to the index by
        # default and we let the normal auto-orientation-code do its thing later
        other_dim = "x" if missing_bar_dim == "y" else "y"
        if not _is_continuous(df_output, args[other_dim]):
            args[missing_bar_dim] = count_name
            df_output = df_output.with_columns(nw.lit(1).alias(count_name))
        else:
            # on the other hand, if the non-missing dimension is continuous, then we
            # can use this information to override the normal auto-orientation code
            if args["orientation"] is None:
                args["orientation"] = "v" if missing_bar_dim == "x" else "h"

    if constructor in hist2d_types:
        del args["orientation"]

    if wide_mode:
        # at this point, `df_output` is semi-long/semi-wide, but we know which columns
        # are which, so we melt it and reassign `args` to refer to the newly-tidy
        # columns, keeping track of various names and manglings set up above
        wide_value_vars = [c for c in args["wide_variable"] if c not in wide_id_vars]
        del args["wide_variable"]
        if wide_cross_name == "__x__":
            wide_cross_name = args["x"]
        elif wide_cross_name == "__y__":
            wide_cross_name = args["y"]
        else:
            wide_cross_name = args["wide_cross"]
        del args["wide_cross"]
        dtype = None
        for v in wide_value_vars:
            v_dtype = df_output.get_column(v).dtype
            v_dtype = "number" if v_dtype.is_numeric() else str(v_dtype)
            if dtype is None:
                dtype = v_dtype
            elif dtype != v_dtype:
                raise ValueError(
                    "Plotly Express cannot process wide-form data with columns of different type."
                )
        df_output = df_output.unpivot(
            index=wide_id_vars,
            on=wide_value_vars,
            variable_name=var_name,
            value_name=value_name,
        )
        assert len(df_output.columns) == len(set(df_output.columns)), (
            "Wide-mode name-inference failure, likely due to a internal bug. "
            "Please report this to "
            "https://github.com/plotly/plotly.py/issues/new and we will try to "
            "replicate and fix it."
        )
        df_output = df_output.with_columns(nw.col(var_name).cast(nw.String))
        orient_v = wide_orientation == "v"

        if hist1d_orientation:
            args["x" if orient_v else "y"] = value_name
            args["y" if orient_v else "x"] = wide_cross_name
            args["color"] = args["color"] or var_name
        elif constructor in [go.Scatter, go.Funnel] + hist2d_types:
            args["x" if orient_v else "y"] = wide_cross_name
            args["y" if orient_v else "x"] = value_name
            if constructor != go.Histogram2d:
                args["color"] = args["color"] or var_name
            if "line_group" in args:
                args["line_group"] = args["line_group"] or var_name
        elif constructor == go.Bar:
            if _is_continuous(df_output, value_name):
                args["x" if orient_v else "y"] = wide_cross_name
                args["y" if orient_v else "x"] = value_name
                args["color"] = args["color"] or var_name
            else:
                args["x" if orient_v else "y"] = value_name
                args["y" if orient_v else "x"] = count_name
                df_output = df_output.with_columns(nw.lit(1).alias(count_name))
                args["color"] = args["color"] or var_name
        elif constructor in [go.Violin, go.Box]:
            args["x" if orient_v else "y"] = wide_cross_name or var_name
            args["y" if orient_v else "x"] = value_name

    if hist1d_orientation and constructor == go.Scatter:
        if args["x"] is not None and args["y"] is not None:
            args["histfunc"] = "sum"
        elif args["x"] is None:
            args["histfunc"] = None
            args["orientation"] = "h"
            args["x"] = count_name
            df_output = df_output.with_columns(nw.lit(1).alias(count_name))
        else:
            args["histfunc"] = None
            args["orientation"] = "v"
            args["y"] = count_name
            df_output = df_output.with_columns(nw.lit(1).alias(count_name))

    if no_color:
        args["color"] = None
    args["data_frame"] = df_output
    return args


def _check_dataframe_all_leaves(df: nw.DataFrame) -> None:
    cols = df.columns
    df_sorted = df.sort(by=cols, descending=False, nulls_last=True)
    null_mask = df_sorted.select(nw.all().is_null())
    df_sorted = df_sorted.select(nw.all().cast(nw.String()))
    null_indices_mask = null_mask.select(
        null_mask=nw.any_horizontal(nw.all())
    ).get_column("null_mask")

    null_mask_filtered = null_mask.filter(null_indices_mask)
    if not null_mask_filtered.is_empty():
        for col_idx in range(1, null_mask_filtered.shape[1]):
            # For each row, if a True value is encountered, then check that
            # all values in subsequent columns are also True
            null_entries_with_non_null_children = (
                ~null_mask_filtered[:, col_idx] & null_mask_filtered[:, col_idx - 1]
            )
            if nw.to_py_scalar(null_entries_with_non_null_children.any()):
                row_idx = null_entries_with_non_null_children.to_list().index(True)
                raise ValueError(
                    "None entries cannot have not-None children",
                    df_sorted.row(row_idx),
                )

    fill_series = nw.new_series(
        name="fill_value",
        values=[""] * len(df_sorted),
        dtype=nw.String(),
        native_namespace=nw.get_native_namespace(df_sorted),
    )
    df_sorted = df_sorted.with_columns(
        **{
            c: df_sorted.get_column(c).zip_with(~null_mask.get_column(c), fill_series)
            for c in cols
        }
    )

    # Conversion to list is due to python native vs pyarrow scalars
    row_strings = (
        df_sorted.select(
            row_strings=nw.concat_str(cols, separator="", ignore_nulls=False)
        )
        .get_column("row_strings")
        .to_list()
    )

    null_indices = set(null_indices_mask.arg_true().to_list())
    for i, (current_row, next_row) in enumerate(
        zip(row_strings[:-1], row_strings[1:]), start=1
    ):
        if (next_row in current_row) and (i in null_indices):
            raise ValueError(
                "Non-leaves rows are not permitted in the dataframe \n",
                df_sorted.row(i),
                "is not a leaf.",
            )


def process_dataframe_hierarchy(args):
    """
    Build dataframe for sunburst, treemap, or icicle when the path argument is provided.
    """
    df: nw.DataFrame = args["data_frame"]
    path = args["path"][::-1]
    _check_dataframe_all_leaves(df[path[::-1]])
    discrete_color = not _is_continuous(df, args["color"]) if args["color"] else False

    df = df.lazy()

    new_path = [col_name + "_path_copy" for col_name in path]
    df = df.with_columns(
        nw.col(col_name).alias(new_col_name)
        for new_col_name, col_name in zip(new_path, path)
    )
    path = new_path
    # ------------ Define aggregation functions --------------------------------
    agg_f = {}
    if args["values"]:
        try:
            df = df.with_columns(nw.col(args["values"]).cast(nw.Float64()))

        except Exception:  # pandas, Polars and pyarrow exception types are different
            raise ValueError(
                "Column `%s` of `df` could not be converted to a numerical data type."
                % args["values"]
            )

        if args["color"] and args["color"] == args["values"]:
            new_value_col_name = args["values"] + "_sum"
            df = df.with_columns(nw.col(args["values"]).alias(new_value_col_name))
            args["values"] = new_value_col_name
        count_colname = args["values"]
    else:
        # we need a count column for the first groupby and the weighted mean of color
        # trick to be sure the col name is unused: take the sum of existing names
        columns = df.collect_schema().names()
        count_colname = (
            "count" if "count" not in columns else "".join([str(el) for el in columns])
        )
        # we can modify df because it's a copy of the px argument
        df = df.with_columns(nw.lit(1).alias(count_colname))
        args["values"] = count_colname

    # Since count_colname is always in agg_f, it can be used later to normalize color
    # in the continuous case after some gymnastic
    agg_f[count_colname] = nw.sum(count_colname)

    discrete_aggs = []
    continuous_aggs = []

    n_unique_token = _generate_temporary_column_name(
        n_bytes=16, columns=df.collect_schema().names()
    )

    # In theory, for discrete columns aggregation, we should have a way to do
    # `.agg(nw.col(x).unique())` in group_by and successively unpack/parse it as:
    # ```
    # (nw.when(nw.col(x).list.len()==1)
    # .then(nw.col(x).list.first())
    # .otherwise(nw.lit("(?)"))
    # )
    # ```
    # which replicates the original pandas only codebase:
    # ```
    # def discrete_agg(x):
    #     uniques = x.unique()
    #     return uniques[0] if len(uniques) == 1 else "(?)"
    #
    # df.groupby(path[i:]).agg(...)
    # ```
    # However this is not possible, therefore the following workaround is provided.
    # We make two aggregations for the same column:
    # - take the max value
    # - take the number of unique values
    # Finally, after the group by statement, it is unpacked via:
    # ```
    # (nw.when(nw.col(col_n_unique) == 1)
    # .then(nw.col(col_max_value))  # which is the unique value
    # .otherwise(nw.lit("(?)"))
    # )
    # ```

    if args["color"]:
        if discrete_color:
            discrete_aggs.append(args["color"])
            agg_f[args["color"]] = nw.col(args["color"]).max()
            agg_f[f"{args['color']}{n_unique_token}"] = (
                nw.col(args["color"])
                .n_unique()
                .alias(f"{args['color']}{n_unique_token}")
            )
        else:
            # This first needs to be multiplied by `count_colname`
            continuous_aggs.append(args["color"])

            agg_f[args["color"]] = nw.sum(args["color"])

    #  Other columns (for color, hover_data, custom_data etc.)
    cols = list(set(df.collect_schema().names()).difference(path))
    df = df.with_columns(nw.col(c).cast(nw.String()) for c in cols if c not in agg_f)

    for col in cols:  # for hover_data, custom_data etc.
        if col not in agg_f:
            # Similar trick as above
            discrete_aggs.append(col)
            agg_f[col] = nw.col(col).max()
            agg_f[f"{col}{n_unique_token}"] = (
                nw.col(col).n_unique().alias(f"{col}{n_unique_token}")
            )
    # Avoid collisions with reserved names - columns in the path have been copied already
    cols = list(set(cols) - set(["labels", "parent", "id"]))
    # ----------------------------------------------------------------------------
    all_trees = []

    if args["color"] and not discrete_color:
        df = df.with_columns(
            (nw.col(args["color"]) * nw.col(count_colname)).alias(args["color"])
        )

    def post_agg(dframe: nw.LazyFrame, continuous_aggs, discrete_aggs) -> nw.LazyFrame:
        """
        - continuous_aggs is either [] or [args["color"]]
        - discrete_aggs is either [args["color"], <rest_of_cols>] or [<rest_of cols>]
        """
        return dframe.with_columns(
            *[nw.col(col) / nw.col(count_colname) for col in continuous_aggs],
            *[
                (
                    nw.when(nw.col(f"{col}{n_unique_token}") == 1)
                    .then(nw.col(col))
                    .otherwise(nw.lit("(?)"))
                    .alias(col)
                )
                for col in discrete_aggs
            ],
        ).drop([f"{col}{n_unique_token}" for col in discrete_aggs])

    for i, level in enumerate(path):
        dfg = (
            df.group_by(path[i:], drop_null_keys=True)
            .agg(**agg_f)
            .pipe(post_agg, continuous_aggs, discrete_aggs)
        )

        # Path label massaging
        df_tree = dfg.with_columns(
            *cols,
            labels=nw.col(level).cast(nw.String()),
            parent=nw.lit(""),
            id=nw.col(level).cast(nw.String()),
        )
        if i < len(path) - 1:
            _concat_str_token = _generate_temporary_column_name(
                n_bytes=16, columns=[*cols, "labels", "parent", "id"]
            )
            df_tree = (
                df_tree.with_columns(
                    nw.concat_str(
                        [
                            nw.col(path[j]).cast(nw.String())
                            for j in range(len(path) - 1, i, -1)
                        ],
                        separator="/",
                    ).alias(_concat_str_token)
                )
                .with_columns(
                    parent=nw.concat_str(
                        [nw.col(_concat_str_token), nw.col("parent")], separator="/"
                    ),
                    id=nw.concat_str(
                        [nw.col(_concat_str_token), nw.col("id")], separator="/"
                    ),
                )
                .drop(_concat_str_token)
            )

        # strip "/" if at the end of the string, equivalent to `.str.rstrip`
        df_tree = df_tree.with_columns(
            parent=nw.col("parent").str.replace("/?$", "").str.replace("^/?", "")
        )

        all_trees.append(df_tree.select(*["labels", "parent", "id", *cols]))

    df_all_trees = nw.maybe_reset_index(nw.concat(all_trees, how="vertical").collect())

    # we want to make sure than (?) is the first color of the sequence
    if args["color"] and discrete_color:
        sort_col_name = "sort_color_if_discrete_color"
        while sort_col_name in df_all_trees.columns:
            sort_col_name += "0"
        df_all_trees = df_all_trees.with_columns(
            nw.col(args["color"]).cast(nw.String()).alias(sort_col_name)
        ).sort(by=sort_col_name, nulls_last=True)

    # Now modify arguments
    args["data_frame"] = df_all_trees
    args["path"] = None
    args["ids"] = "id"
    args["names"] = "labels"
    args["parents"] = "parent"
    if args["color"]:
        if not args["hover_data"]:
            args["hover_data"] = [args["color"]]
        elif isinstance(args["hover_data"], dict):
            if not args["hover_data"].get(args["color"]):
                args["hover_data"][args["color"]] = (True, None)
        else:
            args["hover_data"].append(args["color"])
    return args


def process_dataframe_timeline(args):
    """
    Massage input for bar traces for px.timeline()
    """
    args["is_timeline"] = True
    if args["x_start"] is None or args["x_end"] is None:
        raise ValueError("Both x_start and x_end are required")

    df: nw.DataFrame = args["data_frame"]
    schema = df.schema
    to_convert_to_datetime = [
        col
        for col in [args["x_start"], args["x_end"]]
        if schema[col] != nw.Datetime and schema[col] != nw.Date
    ]

    if to_convert_to_datetime:
        try:
            df = df.with_columns(nw.col(to_convert_to_datetime).str.to_datetime())
        except Exception as exc:
            raise TypeError(
                "Both x_start and x_end must refer to data convertible to datetimes."
            ) from exc

    # note that we are not adding any columns to the data frame here, so no risk of overwrite
    args["data_frame"] = df.with_columns(
        (nw.col(args["x_end"]) - nw.col(args["x_start"]))
        .dt.total_milliseconds()
        .alias(args["x_end"])
    )
    args["x"] = args["x_end"]
    args["base"] = args["x_start"]
    del args["x_start"], args["x_end"]
    return args


def process_dataframe_pie(args, trace_patch):
    import numpy as np

    names = args.get("names")
    if names is None:
        return args, trace_patch
    order_in = args["category_orders"].get(names, {}).copy()
    if not order_in:
        return args, trace_patch
    df: nw.DataFrame = args["data_frame"]
    trace_patch["sort"] = False
    trace_patch["direction"] = "clockwise"
    uniques = df.get_column(names).unique(maintain_order=True).to_list()
    order = [x for x in OrderedDict.fromkeys(list(order_in) + uniques) if x in uniques]

    # Sort args['data_frame'] by column `names` according to order `order`.
    token = nw.generate_temporary_column_name(8, df.columns)
    args["data_frame"] = (
        df.with_columns(
            nw.col(names)
            .replace_strict(order, np.arange(len(order)), return_dtype=nw.UInt32)
            .alias(token)
        )
        .sort(token)
        .drop(token)
    )
    return args, trace_patch


def infer_config(args, constructor, trace_patch, layout_patch):
    attrs = [k for k in direct_attrables + array_attrables if k in args]
    grouped_attrs = []
    df: nw.DataFrame = args["data_frame"]

    # Compute sizeref
    sizeref = 0
    if "size" in args and args["size"]:
        sizeref = (
            nw.to_py_scalar(df.get_column(args["size"]).max()) / args["size_max"] ** 2
        )

    # Compute color attributes and grouping attributes
    if "color" in args:
        if "color_continuous_scale" in args:
            if "color_discrete_sequence" not in args:
                attrs.append("color")
            else:
                if args["color"] and _is_continuous(df, args["color"]):
                    attrs.append("color")
                    args["color_is_continuous"] = True
                elif constructor in [go.Sunburst, go.Treemap, go.Icicle]:
                    attrs.append("color")
                    args["color_is_continuous"] = False
                else:
                    grouped_attrs.append("marker.color")
        elif "line_group" in args or constructor == go.Histogram2dContour:
            grouped_attrs.append("line.color")
        elif constructor in [go.Pie, go.Funnelarea]:
            attrs.append("color")
            if args["color"]:
                if args["hover_data"] is None:
                    args["hover_data"] = []
                args["hover_data"].append(args["color"])
        else:
            grouped_attrs.append("marker.color")

        show_colorbar = bool(
            "color" in attrs
            and args["color"]
            and constructor not in [go.Pie, go.Funnelarea]
            and (
                constructor not in [go.Treemap, go.Sunburst, go.Icicle]
                or args.get("color_is_continuous")
            )
        )
    else:
        show_colorbar = False

    if "line_dash" in args:
        grouped_attrs.append("line.dash")

    if "symbol" in args:
        grouped_attrs.append("marker.symbol")

    if "pattern_shape" in args:
        if constructor in [go.Scatter]:
            grouped_attrs.append("fillpattern.shape")
        else:
            grouped_attrs.append("marker.pattern.shape")

    if "orientation" in args:
        has_x = args["x"] is not None
        has_y = args["y"] is not None
        if args["orientation"] is None:
            if constructor in [go.Histogram, go.Scatter]:
                if has_y and not has_x:
                    args["orientation"] = "h"
            elif constructor in [go.Violin, go.Box, go.Bar, go.Funnel]:
                if has_x and not has_y:
                    args["orientation"] = "h"

        if args["orientation"] is None and has_x and has_y:
            x_is_continuous = _is_continuous(df, args["x"])
            y_is_continuous = _is_continuous(df, args["y"])
            if x_is_continuous and not y_is_continuous:
                args["orientation"] = "h"
            if y_is_continuous and not x_is_continuous:
                args["orientation"] = "v"

        if args["orientation"] is None:
            args["orientation"] = "v"

        if constructor == go.Histogram:
            if has_x and has_y and args["histfunc"] is None:
                args["histfunc"] = trace_patch["histfunc"] = "sum"

            orientation = args["orientation"]
            nbins = args["nbins"]
            trace_patch["nbinsx"] = nbins if orientation == "v" else None
            trace_patch["nbinsy"] = None if orientation == "v" else nbins
            trace_patch["bingroup"] = "x" if orientation == "v" else "y"
        trace_patch["orientation"] = args["orientation"]

        if constructor in [go.Violin, go.Box]:
            mode = "boxmode" if constructor == go.Box else "violinmode"
            if layout_patch[mode] is None and args["color"] is not None:
                if args["y"] == args["color"] and args["orientation"] == "h":
                    layout_patch[mode] = "overlay"
                elif args["x"] == args["color"] and args["orientation"] == "v":
                    layout_patch[mode] = "overlay"
            if layout_patch[mode] is None:
                layout_patch[mode] = "group"

    if (
        constructor == go.Histogram2d
        and args["z"] is not None
        and args["histfunc"] is None
    ):
        args["histfunc"] = trace_patch["histfunc"] = "sum"

    if args.get("text_auto", False) is not False:
        if constructor in [go.Histogram2d, go.Histogram2dContour]:
            letter = "z"
        elif constructor == go.Bar:
            letter = "y" if args["orientation"] == "v" else "x"
        else:
            letter = "value"
        if args["text_auto"] is True:
            trace_patch["texttemplate"] = "%{" + letter + "}"
        else:
            trace_patch["texttemplate"] = "%{" + letter + ":" + args["text_auto"] + "}"

    if constructor in [go.Histogram2d, go.Densitymap, go.Densitymapbox]:
        show_colorbar = True
        trace_patch["coloraxis"] = "coloraxis1"

    if "opacity" in args:
        if args["opacity"] is None:
            if "barmode" in args and args["barmode"] == "overlay":
                trace_patch["marker"] = dict(opacity=0.5)
        elif constructor in [
            go.Densitymap,
            go.Densitymapbox,
            go.Pie,
            go.Funnel,
            go.Funnelarea,
        ]:
            trace_patch["opacity"] = args["opacity"]
        else:
            trace_patch["marker"] = dict(opacity=args["opacity"])
    if (
        "line_group" in args or "line_dash" in args
    ):  # px.line, px.line_*, px.area, px.ecdf
        modes = set()
        if args.get("lines", True):
            modes.add("lines")
        if args.get("text") or args.get("symbol") or args.get("markers"):
            modes.add("markers")
        if args.get("text"):
            modes.add("text")
        if len(modes) == 0:
            modes.add("lines")
        trace_patch["mode"] = "+".join(sorted(modes))
    elif constructor != go.Splom and (
        "symbol" in args or constructor in [go.Scattermap, go.Scattermapbox]
    ):
        trace_patch["mode"] = "markers" + ("+text" if args["text"] else "")

    if "line_shape" in args:
        trace_patch["line"] = dict(shape=args["line_shape"])
    elif "ecdfmode" in args:
        trace_patch["line"] = dict(
            shape="vh" if args["ecdfmode"] == "reversed" else "hv"
        )

    if "geojson" in args:
        trace_patch["featureidkey"] = args["featureidkey"]
        trace_patch["geojson"] = (
            args["geojson"]
            if not hasattr(args["geojson"], "__geo_interface__")  # for geopandas
            else args["geojson"].__geo_interface__
        )

    # Compute marginal attribute: copy to appropriate marginal_*
    if "marginal" in args:
        position = "marginal_x" if args["orientation"] == "v" else "marginal_y"
        other_position = "marginal_x" if args["orientation"] == "h" else "marginal_y"
        args[position] = args["marginal"]
        args[other_position] = None

    # Ignore facet rows and columns when data frame is empty so as to prevent nrows/ncols equaling 0
    if df.is_empty():
        args["facet_row"] = args["facet_col"] = None

    # If both marginals and faceting are specified, faceting wins
    if args.get("facet_col") is not None and args.get("marginal_y") is not None:
        args["marginal_y"] = None

    if args.get("facet_row") is not None and args.get("marginal_x") is not None:
        args["marginal_x"] = None

    # facet_col_wrap only works if no marginals or row faceting is used
    if (
        args.get("marginal_x") is not None
        or args.get("marginal_y") is not None
        or args.get("facet_row") is not None
    ):
        args["facet_col_wrap"] = 0

    if "trendline" in args and args["trendline"] is not None:
        if args["trendline"] not in trendline_functions:
            raise ValueError(
                "Value '%s' for `trendline` must be one of %s"
                % (args["trendline"], trendline_functions.keys())
            )

    if "trendline_options" in args and args["trendline_options"] is None:
        args["trendline_options"] = dict()

    if "ecdfnorm" in args:
        if args.get("ecdfnorm", None) not in [None, "percent", "probability"]:
            raise ValueError(
                "`ecdfnorm` must be one of None, 'percent' or 'probability'. "
                + "'%s' was provided." % args["ecdfnorm"]
            )
        args["histnorm"] = args["ecdfnorm"]

    # Compute applicable grouping attributes
    grouped_attrs.extend([k for k in group_attrables if k in args])

    # Create grouped mappings
    grouped_mappings = [make_mapping(args, a) for a in grouped_attrs]

    # Create trace specs
    trace_specs = make_trace_spec(args, constructor, attrs, trace_patch)
    return trace_specs, grouped_mappings, sizeref, show_colorbar


def get_groups_and_orders(args, grouper):
    """
    `orders` is the user-supplied ordering with the remaining data-frame-supplied
    ordering appended if the column is used for grouping. It includes anything the user
    gave, for any variable, including values not present in the dataset. It's a dict
    where the keys are e.g. "x" or "color"

    `groups` is the dicts of groups, ordered by the order above. Its keys are
    tuples like [("value1", ""), ("value2", "")] where each tuple contains the name
    of a single dimension-group
    """
    orders = {} if "category_orders" not in args else args["category_orders"].copy()
    df: nw.DataFrame = args["data_frame"]
    # figure out orders and what the single group name would be if there were one
    single_group_name = []
    unique_cache = dict()

    for i, col in enumerate(grouper):
        if col == one_group:
            single_group_name.append("")
        else:
            if col not in unique_cache:
                unique_cache[col] = (
                    df.get_column(col).unique(maintain_order=True).to_list()
                )
            uniques = unique_cache[col]
            if len(uniques) == 1:
                single_group_name.append(uniques[0])
            if col not in orders:
                orders[col] = uniques
            else:
                orders[col] = list(OrderedDict.fromkeys(list(orders[col]) + uniques))

    if len(single_group_name) == len(grouper):
        # we have a single group, so we can skip all group-by operations!
        groups = {tuple(single_group_name): df}
    else:
        required_grouper = [group for group in orders if group in grouper]
        grouped = dict(df.group_by(required_grouper, drop_null_keys=True).__iter__())

        sorted_group_names = sorted(
            grouped.keys(),
            key=lambda values: [
                orders[group].index(value) if value in orders[group] else -1
                for group, value in zip(required_grouper, values)
            ],
        )

        # calculate the full group_names by inserting "" in the tuple index for one_group groups
        full_sorted_group_names = [
            tuple(
                [
                    (
                        ""
                        if col == one_group
                        else sub_group_names[required_grouper.index(col)]
                    )
                    for col in grouper
                ]
            )
            for sub_group_names in sorted_group_names
        ]

        groups = {
            sf: grouped[s] for sf, s in zip(full_sorted_group_names, sorted_group_names)
        }
    return groups, orders


def make_figure(args, constructor, trace_patch=None, layout_patch=None):
    trace_patch = trace_patch or {}
    layout_patch = layout_patch or {}
    apply_default_cascade(args)

    args = build_dataframe(args, constructor)
    if constructor in [go.Treemap, go.Sunburst, go.Icicle] and args["path"] is not None:
        args = process_dataframe_hierarchy(args)
    if constructor in [go.Pie]:
        args, trace_patch = process_dataframe_pie(args, trace_patch)
    if constructor == "timeline":
        constructor = go.Bar
        args = process_dataframe_timeline(args)

    # If we have marginal histograms, set barmode to "overlay"
    if "histogram" in [args.get("marginal_x"), args.get("marginal_y")]:
        layout_patch["barmode"] = "overlay"

    trace_specs, grouped_mappings, sizeref, show_colorbar = infer_config(
        args, constructor, trace_patch, layout_patch
    )
    grouper = [x.grouper or one_group for x in grouped_mappings] or [one_group]
    groups, orders = get_groups_and_orders(args, grouper)

    col_labels = []
    row_labels = []
    nrows = ncols = 1
    for m in grouped_mappings:
        if m.grouper not in orders:
            m.val_map[""] = m.sequence[0]
        else:
            sorted_values = orders[m.grouper]
            if m.facet == "col":
                prefix = get_label(args, args["facet_col"]) + "="
                col_labels = [prefix + str(s) for s in sorted_values]
                ncols = len(col_labels)
            if m.facet == "row":
                prefix = get_label(args, args["facet_row"]) + "="
                row_labels = [prefix + str(s) for s in sorted_values]
                nrows = len(row_labels)
            for val in sorted_values:
                if val not in m.val_map:  # always False if it's an IdentityMap
                    m.val_map[val] = m.sequence[len(m.val_map) % len(m.sequence)]

    subplot_type = _subplot_type_for_trace_type(constructor().type)

    trace_names_by_frame = {}
    frames = OrderedDict()
    trendline_rows = []
    trace_name_labels = None
    facet_col_wrap = args.get("facet_col_wrap", 0)
    for group_name, group in groups.items():
        mapping_labels = OrderedDict()
        trace_name_labels = OrderedDict()
        frame_name = ""
        for col, val, m in zip(grouper, group_name, grouped_mappings):
            if col != one_group:
                key = get_label(args, col)
                if not isinstance(m.val_map, IdentityMap):
                    mapping_labels[key] = str(val)
                    if m.show_in_trace_name:
                        trace_name_labels[key] = str(val)
                if m.variable == "animation_frame":
                    frame_name = val
        trace_name = ", ".join(trace_name_labels.values())
        if frame_name not in trace_names_by_frame:
            trace_names_by_frame[frame_name] = set()
        trace_names = trace_names_by_frame[frame_name]

        for trace_spec in trace_specs:
            # Create the trace
            trace = trace_spec.constructor(name=trace_name)
            if trace_spec.constructor not in [
                go.Parcats,
                go.Parcoords,
                go.Choropleth,
                go.Choroplethmap,
                go.Choroplethmapbox,
                go.Densitymap,
                go.Densitymapbox,
                go.Histogram2d,
                go.Sunburst,
                go.Treemap,
                go.Icicle,
            ]:
                trace.update(
                    legendgroup=trace_name,
                    showlegend=(trace_name != "" and trace_name not in trace_names),
                )

            # Set 'offsetgroup' only in group barmode (or if no barmode is set)
            barmode = layout_patch.get("barmode")
            if trace_spec.constructor in [go.Bar, go.Box, go.Violin, go.Histogram] and (
                barmode == "group" or barmode is None
            ):
                trace.update(alignmentgroup=True, offsetgroup=trace_name)
            trace_names.add(trace_name)

            # Init subplot row/col
            trace._subplot_row = 1
            trace._subplot_col = 1

            for i, m in enumerate(grouped_mappings):
                val = group_name[i]
                try:
                    m.updater(trace, m.val_map[val])  # covers most cases
                except ValueError:
                    # this catches some odd cases like marginals
                    if (
                        trace_spec != trace_specs[0]
                        and (
                            trace_spec.constructor in [go.Violin, go.Box]
                            and m.variable in ["symbol", "pattern", "dash"]
                        )
                        or (
                            trace_spec.constructor in [go.Histogram]
                            and m.variable in ["symbol", "dash"]
                        )
                    ):
                        pass
                    elif (
                        trace_spec != trace_specs[0]
                        and trace_spec.constructor in [go.Histogram]
                        and m.variable == "color"
                    ):
                        trace.update(marker=dict(color=m.val_map[val]))
                    elif (
                        trace_spec.constructor
                        in [go.Choropleth, go.Choroplethmap, go.Choroplethmapbox]
                        and m.variable == "color"
                    ):
                        trace.update(
                            z=[1] * len(group),
                            colorscale=[m.val_map[val]] * 2,
                            showscale=False,
                            showlegend=True,
                        )
                    else:
                        raise

                # Find row for trace, handling facet_row and marginal_x
                if m.facet == "row":
                    row = m.val_map[val]
                else:
                    if (
                        args.get("marginal_x") is not None  # there is a marginal
                        and trace_spec.marginal != "x"  # and we're not it
                    ):
                        row = 2
                    else:
                        row = 1

                # Find col for trace, handling facet_col and marginal_y
                if m.facet == "col":
                    col = m.val_map[val]
                    if facet_col_wrap:  # assumes no facet_row, no marginals
                        row = 1 + ((col - 1) // facet_col_wrap)
                        col = 1 + ((col - 1) % facet_col_wrap)
                else:
                    if trace_spec.marginal == "y":
                        col = 2
                    else:
                        col = 1

                if row > 1:
                    trace._subplot_row = row

                if col > 1:
                    trace._subplot_col = col
            if (
                trace_specs[0].constructor == go.Histogram2dContour
                and trace_spec.constructor == go.Box
                and trace.line.color
            ):
                trace.update(marker=dict(color=trace.line.color))

            if "ecdfmode" in args:
                base = args["x"] if args["orientation"] == "v" else args["y"]
                var = args["x"] if args["orientation"] == "h" else args["y"]
                ascending = args.get("ecdfmode", "standard") != "reversed"
                group = group.sort(by=base, descending=not ascending, nulls_last=True)
                group_sum = group.get_column(
                    var
                ).sum()  # compute here before next line mutates
                group = group.with_columns(nw.col(var).cum_sum().alias(var))
                if not ascending:
                    group = group.sort(by=base, descending=False, nulls_last=True)

                if args.get("ecdfmode", "standard") == "complementary":
                    group = group.with_columns((group_sum - nw.col(var)).alias(var))

                if args["ecdfnorm"] == "probability":
                    group = group.with_columns(nw.col(var) / group_sum)
                elif args["ecdfnorm"] == "percent":
                    group = group.with_columns((nw.col(var) / group_sum) * 100.0)

            patch, fit_results = make_trace_kwargs(
                args, trace_spec, group, mapping_labels.copy(), sizeref
            )
            trace.update(patch)
            if fit_results is not None:
                trendline_rows.append(mapping_labels.copy())
                trendline_rows[-1]["px_fit_results"] = fit_results
            if frame_name not in frames:
                frames[frame_name] = dict(data=[], name=frame_name)
            frames[frame_name]["data"].append(trace)
    frame_list = [f for f in frames.values()]
    if len(frame_list) > 1:
        frame_list = sorted(
            frame_list, key=lambda f: orders[args["animation_frame"]].index(f["name"])
        )

    if show_colorbar:
        colorvar = (
            "z"
            if constructor in [go.Histogram2d, go.Densitymap, go.Densitymapbox]
            else "color"
        )
        range_color = args["range_color"] or [None, None]

        colorscale_validator = ColorscaleValidator("colorscale", "make_figure")
        layout_patch["coloraxis1"] = dict(
            colorscale=colorscale_validator.validate_coerce(
                args["color_continuous_scale"]
            ),
            cmid=args["color_continuous_midpoint"],
            cmin=range_color[0],
            cmax=range_color[1],
            colorbar=dict(
                title_text=get_decorated_label(args, args[colorvar], colorvar)
            ),
        )
    for v in ["height", "width"]:
        if args[v]:
            layout_patch[v] = args[v]
    layout_patch["legend"] = dict(tracegroupgap=0)
    if trace_name_labels:
        layout_patch["legend"]["title_text"] = ", ".join(trace_name_labels)
    if args["title"]:
        layout_patch["title_text"] = args["title"]
    elif args["template"].layout.margin.t is None:
        layout_patch["margin"] = {"t": 60}
    if args["subtitle"]:
        layout_patch["title_subtitle_text"] = args["subtitle"]
    if (
        "size" in args
        and args["size"]
        and args["template"].layout.legend.itemsizing is None
    ):
        layout_patch["legend"]["itemsizing"] = "constant"

    if facet_col_wrap:
        nrows = math.ceil(ncols / facet_col_wrap)
        ncols = min(ncols, facet_col_wrap)

    if args.get("marginal_x") is not None:
        nrows += 1

    if args.get("marginal_y") is not None:
        ncols += 1

    fig = init_figure(
        args, subplot_type, frame_list, nrows, ncols, col_labels, row_labels
    )

    # Position traces in subplots
    for frame in frame_list:
        for trace in frame["data"]:
            if isinstance(trace, go.Splom):
                # Special case that is not compatible with make_subplots
                continue

            _set_trace_grid_reference(
                trace,
                fig.layout,
                fig._grid_ref,
                nrows - trace._subplot_row + 1,
                trace._subplot_col,
            )

    # Add traces, layout and frames to figure
    fig.add_traces(frame_list[0]["data"] if len(frame_list) > 0 else [])
    fig.update_layout(layout_patch)
    if "template" in args and args["template"] is not None:
        fig.update_layout(template=args["template"], overwrite=True)
    for f in frame_list:
        f["name"] = str(f["name"])
    fig.frames = frame_list if len(frames) > 1 else []

    if args.get("trendline") and args.get("trendline_scope", "trace") == "overall":
        trendline_spec = make_trendline_spec(args, constructor)
        trendline_trace = trendline_spec.constructor(
            name="Overall Trendline", legendgroup="Overall Trendline", showlegend=False
        )
        if "line" not in trendline_spec.trace_patch:  # no color override
            for m in grouped_mappings:
                if m.variable == "color":
                    next_color = m.sequence[len(m.val_map) % len(m.sequence)]
                    trendline_spec.trace_patch["line"] = dict(color=next_color)
        patch, fit_results = make_trace_kwargs(
            args, trendline_spec, args["data_frame"], {}, sizeref
        )
        trendline_trace.update(patch)
        fig.add_trace(
            trendline_trace, row="all", col="all", exclude_empty_subplots=True
        )
        fig.update_traces(selector=-1, showlegend=True)
        if fit_results is not None:
            trendline_rows.append(dict(px_fit_results=fit_results))

    if trendline_rows:
        try:
            import pandas as pd

            fig._px_trendlines = pd.DataFrame(trendline_rows)
        except ImportError:
            msg = "Trendlines require pandas to be installed."
            raise NotImplementedError(msg)
    else:
        fig._px_trendlines = []

    configure_axes(args, constructor, fig, orders)
    configure_animation_controls(args, constructor, fig)
    return fig


def init_figure(args, subplot_type, frame_list, nrows, ncols, col_labels, row_labels):
    # Build subplot specs
    specs = [[dict(type=subplot_type or "domain")] * ncols for _ in range(nrows)]

    # Default row/column widths uniform
    column_widths = [1.0] * ncols
    row_heights = [1.0] * nrows
    facet_col_wrap = args.get("facet_col_wrap", 0)

    # Build column_widths/row_heights
    if subplot_type == "xy":
        if args.get("marginal_x") is not None:
            if args["marginal_x"] == "histogram" or ("color" in args and args["color"]):
                main_size = 0.74
            else:
                main_size = 0.84

            row_heights = [main_size] * (nrows - 1) + [1 - main_size]
            vertical_spacing = 0.01
        elif facet_col_wrap:
            vertical_spacing = args.get("facet_row_spacing") or 0.07
        else:
            vertical_spacing = args.get("facet_row_spacing") or 0.03

        if args.get("marginal_y") is not None:
            if args["marginal_y"] == "histogram" or ("color" in args and args["color"]):
                main_size = 0.74
            else:
                main_size = 0.84

            column_widths = [main_size] * (ncols - 1) + [1 - main_size]
            horizontal_spacing = 0.005
        else:
            horizontal_spacing = args.get("facet_col_spacing") or 0.02
    else:
        # Other subplot types:
        #   'scene', 'geo', 'polar', 'ternary', 'mapbox', 'domain', None
        #
        # We can customize subplot spacing per type once we enable faceting
        # for all plot types
        if facet_col_wrap:
            vertical_spacing = args.get("facet_row_spacing") or 0.07
        else:
            vertical_spacing = args.get("facet_row_spacing") or 0.03
        horizontal_spacing = args.get("facet_col_spacing") or 0.02

    if facet_col_wrap:
        subplot_labels = [None] * nrows * ncols
        while len(col_labels) < nrows * ncols:
            col_labels.append(None)
        for i in range(nrows):
            for j in range(ncols):
                subplot_labels[i * ncols + j] = col_labels[(nrows - 1 - i) * ncols + j]

    def _spacing_error_translator(e, direction, facet_arg):
        """
        Translates the spacing errors thrown by the underlying make_subplots
        routine into one that describes an argument adjustable through px.
        """
        if ("%s spacing" % (direction,)) in e.args[0]:
            e.args = (
                e.args[0]
                + """
Use the {facet_arg} argument to adjust this spacing.""".format(facet_arg=facet_arg),
            )
            raise e

    # Create figure with subplots
    try:
        fig = make_subplots(
            rows=nrows,
            cols=ncols,
            specs=specs,
            shared_xaxes="all",
            shared_yaxes="all",
            row_titles=[] if facet_col_wrap else list(reversed(row_labels)),
            column_titles=[] if facet_col_wrap else col_labels,
            subplot_titles=subplot_labels if facet_col_wrap else [],
            horizontal_spacing=horizontal_spacing,
            vertical_spacing=vertical_spacing,
            row_heights=row_heights,
            column_widths=column_widths,
            start_cell="bottom-left",
        )
    except ValueError as e:
        _spacing_error_translator(e, "Horizontal", "facet_col_spacing")
        _spacing_error_translator(e, "Vertical", "facet_row_spacing")
        raise

    # Remove explicit font size of row/col titles so template can take over
    for annot in fig.layout.annotations:
        annot.update(font=None)

    return fig
