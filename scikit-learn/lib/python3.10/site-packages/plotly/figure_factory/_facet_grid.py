from plotly import exceptions, optional_imports
import plotly.colors as clrs
from plotly.figure_factory import utils
from plotly.subplots import make_subplots

import math
from numbers import Number

pd = optional_imports.get_module("pandas")

TICK_COLOR = "#969696"
AXIS_TITLE_COLOR = "#0f0f0f"
AXIS_TITLE_SIZE = 12
GRID_COLOR = "#ffffff"
LEGEND_COLOR = "#efefef"
PLOT_BGCOLOR = "#ededed"
ANNOT_RECT_COLOR = "#d0d0d0"
LEGEND_BORDER_WIDTH = 1
LEGEND_ANNOT_X = 1.05
LEGEND_ANNOT_Y = 0.5
MAX_TICKS_PER_AXIS = 5
THRES_FOR_FLIPPED_FACET_TITLES = 10
GRID_WIDTH = 1

VALID_TRACE_TYPES = ["scatter", "scattergl", "histogram", "bar", "box"]

CUSTOM_LABEL_ERROR = (
    "If you are using a dictionary for custom labels for the facet row/col, "
    "make sure each key in that column of the dataframe is in your facet "
    "labels. The keys you need are {}"
)


def _is_flipped(num):
    if num >= THRES_FOR_FLIPPED_FACET_TITLES:
        flipped = True
    else:
        flipped = False
    return flipped


def _return_label(original_label, facet_labels, facet_var):
    if isinstance(facet_labels, dict):
        label = facet_labels[original_label]
    elif isinstance(facet_labels, str):
        label = "{}: {}".format(facet_var, original_label)
    else:
        label = original_label
    return label


def _legend_annotation(color_name):
    legend_title = dict(
        textangle=0,
        xanchor="left",
        yanchor="middle",
        x=LEGEND_ANNOT_X,
        y=1.03,
        showarrow=False,
        xref="paper",
        yref="paper",
        text="factor({})".format(color_name),
        font=dict(size=13, color="#000000"),
    )
    return legend_title


def _annotation_dict(
    text, lane, num_of_lanes, SUBPLOT_SPACING, row_col="col", flipped=True
):
    temp = (1 - (num_of_lanes - 1) * SUBPLOT_SPACING) / (num_of_lanes)
    if not flipped:
        xanchor = "center"
        yanchor = "middle"
        if row_col == "col":
            x = (lane - 1) * (temp + SUBPLOT_SPACING) + 0.5 * temp
            y = 1.03
            textangle = 0
        elif row_col == "row":
            y = (lane - 1) * (temp + SUBPLOT_SPACING) + 0.5 * temp
            x = 1.03
            textangle = 90
    else:
        if row_col == "col":
            xanchor = "center"
            yanchor = "bottom"
            x = (lane - 1) * (temp + SUBPLOT_SPACING) + 0.5 * temp
            y = 1.0
            textangle = 270
        elif row_col == "row":
            xanchor = "left"
            yanchor = "middle"
            y = (lane - 1) * (temp + SUBPLOT_SPACING) + 0.5 * temp
            x = 1.0
            textangle = 0

    annotation_dict = dict(
        textangle=textangle,
        xanchor=xanchor,
        yanchor=yanchor,
        x=x,
        y=y,
        showarrow=False,
        xref="paper",
        yref="paper",
        text=str(text),
        font=dict(size=13, color=AXIS_TITLE_COLOR),
    )
    return annotation_dict


def _axis_title_annotation(text, x_or_y_axis):
    if x_or_y_axis == "x":
        x_pos = 0.5
        y_pos = -0.1
        textangle = 0
    elif x_or_y_axis == "y":
        x_pos = -0.1
        y_pos = 0.5
        textangle = 270

    if not text:
        text = ""

    annot = {
        "font": {"color": "#000000", "size": AXIS_TITLE_SIZE},
        "showarrow": False,
        "text": text,
        "textangle": textangle,
        "x": x_pos,
        "xanchor": "center",
        "xref": "paper",
        "y": y_pos,
        "yanchor": "middle",
        "yref": "paper",
    }
    return annot


def _add_shapes_to_fig(fig, annot_rect_color, flipped_rows=False, flipped_cols=False):
    shapes_list = []
    for key in fig["layout"].to_plotly_json().keys():
        if "axis" in key and fig["layout"][key]["domain"] != [0.0, 1.0]:
            shape = {
                "fillcolor": annot_rect_color,
                "layer": "below",
                "line": {"color": annot_rect_color, "width": 1},
                "type": "rect",
                "xref": "paper",
                "yref": "paper",
            }

            if "xaxis" in key:
                shape["x0"] = fig["layout"][key]["domain"][0]
                shape["x1"] = fig["layout"][key]["domain"][1]
                shape["y0"] = 1.005
                shape["y1"] = 1.05

                if flipped_cols:
                    shape["y1"] += 0.5
                shapes_list.append(shape)

            elif "yaxis" in key:
                shape["x0"] = 1.005
                shape["x1"] = 1.05
                shape["y0"] = fig["layout"][key]["domain"][0]
                shape["y1"] = fig["layout"][key]["domain"][1]

                if flipped_rows:
                    shape["x1"] += 1
                shapes_list.append(shape)

    fig["layout"]["shapes"] = shapes_list


def _make_trace_for_scatter(trace, trace_type, color, **kwargs_marker):
    if trace_type in ["scatter", "scattergl"]:
        trace["mode"] = "markers"
        trace["marker"] = dict(color=color, **kwargs_marker)
    return trace


def _facet_grid_color_categorical(
    df,
    x,
    y,
    facet_row,
    facet_col,
    color_name,
    colormap,
    num_of_rows,
    num_of_cols,
    facet_row_labels,
    facet_col_labels,
    trace_type,
    flipped_rows,
    flipped_cols,
    show_boxes,
    SUBPLOT_SPACING,
    marker_color,
    kwargs_trace,
    kwargs_marker,
):
    fig = make_subplots(
        rows=num_of_rows,
        cols=num_of_cols,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=SUBPLOT_SPACING,
        vertical_spacing=SUBPLOT_SPACING,
        print_grid=False,
    )

    annotations = []
    if not facet_row and not facet_col:
        color_groups = list(df.groupby(color_name))
        for group in color_groups:
            trace = dict(
                type=trace_type,
                name=group[0],
                marker=dict(color=colormap[group[0]]),
                **kwargs_trace,
            )
            if x:
                trace["x"] = group[1][x]
            if y:
                trace["y"] = group[1][y]
            trace = _make_trace_for_scatter(
                trace, trace_type, colormap[group[0]], **kwargs_marker
            )

            fig.append_trace(trace, 1, 1)

    elif (facet_row and not facet_col) or (not facet_row and facet_col):
        groups_by_facet = list(df.groupby(facet_row if facet_row else facet_col))
        for j, group in enumerate(groups_by_facet):
            for color_val in df[color_name].unique():
                data_by_color = group[1][group[1][color_name] == color_val]
                trace = dict(
                    type=trace_type,
                    name=color_val,
                    marker=dict(color=colormap[color_val]),
                    **kwargs_trace,
                )
                if x:
                    trace["x"] = data_by_color[x]
                if y:
                    trace["y"] = data_by_color[y]
                trace = _make_trace_for_scatter(
                    trace, trace_type, colormap[color_val], **kwargs_marker
                )

                fig.append_trace(
                    trace, j + 1 if facet_row else 1, 1 if facet_row else j + 1
                )

            label = _return_label(
                group[0],
                facet_row_labels if facet_row else facet_col_labels,
                facet_row if facet_row else facet_col,
            )

            annotations.append(
                _annotation_dict(
                    label,
                    num_of_rows - j if facet_row else j + 1,
                    num_of_rows if facet_row else num_of_cols,
                    SUBPLOT_SPACING,
                    "row" if facet_row else "col",
                    flipped_rows,
                )
            )

    elif facet_row and facet_col:
        groups_by_facets = list(df.groupby([facet_row, facet_col]))
        tuple_to_facet_group = {item[0]: item[1] for item in groups_by_facets}

        row_values = df[facet_row].unique()
        col_values = df[facet_col].unique()
        color_vals = df[color_name].unique()
        for row_count, x_val in enumerate(row_values):
            for col_count, y_val in enumerate(col_values):
                try:
                    group = tuple_to_facet_group[(x_val, y_val)]
                except KeyError:
                    group = pd.DataFrame(
                        [[None, None, None]], columns=[x, y, color_name]
                    )

                for color_val in color_vals:
                    if group.values.tolist() != [[None, None, None]]:
                        group_filtered = group[group[color_name] == color_val]

                        trace = dict(
                            type=trace_type,
                            name=color_val,
                            marker=dict(color=colormap[color_val]),
                            **kwargs_trace,
                        )
                        new_x = group_filtered[x]
                        new_y = group_filtered[y]
                    else:
                        trace = dict(
                            type=trace_type,
                            name=color_val,
                            marker=dict(color=colormap[color_val]),
                            showlegend=False,
                            **kwargs_trace,
                        )
                        new_x = group[x]
                        new_y = group[y]

                    if x:
                        trace["x"] = new_x
                    if y:
                        trace["y"] = new_y
                    trace = _make_trace_for_scatter(
                        trace, trace_type, colormap[color_val], **kwargs_marker
                    )

                    fig.append_trace(trace, row_count + 1, col_count + 1)
                if row_count == 0:
                    label = _return_label(
                        col_values[col_count], facet_col_labels, facet_col
                    )
                    annotations.append(
                        _annotation_dict(
                            label,
                            col_count + 1,
                            num_of_cols,
                            SUBPLOT_SPACING,
                            row_col="col",
                            flipped=flipped_cols,
                        )
                    )
            label = _return_label(row_values[row_count], facet_row_labels, facet_row)
            annotations.append(
                _annotation_dict(
                    label,
                    num_of_rows - row_count,
                    num_of_rows,
                    SUBPLOT_SPACING,
                    row_col="row",
                    flipped=flipped_rows,
                )
            )

    return fig, annotations


def _facet_grid_color_numerical(
    df,
    x,
    y,
    facet_row,
    facet_col,
    color_name,
    colormap,
    num_of_rows,
    num_of_cols,
    facet_row_labels,
    facet_col_labels,
    trace_type,
    flipped_rows,
    flipped_cols,
    show_boxes,
    SUBPLOT_SPACING,
    marker_color,
    kwargs_trace,
    kwargs_marker,
):
    fig = make_subplots(
        rows=num_of_rows,
        cols=num_of_cols,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=SUBPLOT_SPACING,
        vertical_spacing=SUBPLOT_SPACING,
        print_grid=False,
    )

    annotations = []
    if not facet_row and not facet_col:
        trace = dict(
            type=trace_type,
            marker=dict(color=df[color_name], colorscale=colormap, showscale=True),
            **kwargs_trace,
        )
        if x:
            trace["x"] = df[x]
        if y:
            trace["y"] = df[y]
        trace = _make_trace_for_scatter(
            trace, trace_type, df[color_name], **kwargs_marker
        )

        fig.append_trace(trace, 1, 1)

    if (facet_row and not facet_col) or (not facet_row and facet_col):
        groups_by_facet = list(df.groupby(facet_row if facet_row else facet_col))
        for j, group in enumerate(groups_by_facet):
            trace = dict(
                type=trace_type,
                marker=dict(
                    color=df[color_name],
                    colorscale=colormap,
                    showscale=True,
                    colorbar=dict(x=1.15),
                ),
                **kwargs_trace,
            )
            if x:
                trace["x"] = group[1][x]
            if y:
                trace["y"] = group[1][y]
            trace = _make_trace_for_scatter(
                trace, trace_type, df[color_name], **kwargs_marker
            )

            fig.append_trace(
                trace, j + 1 if facet_row else 1, 1 if facet_row else j + 1
            )

            labels = facet_row_labels if facet_row else facet_col_labels
            label = _return_label(
                group[0], labels, facet_row if facet_row else facet_col
            )

            annotations.append(
                _annotation_dict(
                    label,
                    num_of_rows - j if facet_row else j + 1,
                    num_of_rows if facet_row else num_of_cols,
                    SUBPLOT_SPACING,
                    "row" if facet_row else "col",
                    flipped=flipped_rows,
                )
            )

    elif facet_row and facet_col:
        groups_by_facets = list(df.groupby([facet_row, facet_col]))
        tuple_to_facet_group = {item[0]: item[1] for item in groups_by_facets}

        row_values = df[facet_row].unique()
        col_values = df[facet_col].unique()
        for row_count, x_val in enumerate(row_values):
            for col_count, y_val in enumerate(col_values):
                try:
                    group = tuple_to_facet_group[(x_val, y_val)]
                except KeyError:
                    group = pd.DataFrame(
                        [[None, None, None]], columns=[x, y, color_name]
                    )

                if group.values.tolist() != [[None, None, None]]:
                    trace = dict(
                        type=trace_type,
                        marker=dict(
                            color=df[color_name],
                            colorscale=colormap,
                            showscale=(row_count == 0),
                            colorbar=dict(x=1.15),
                        ),
                        **kwargs_trace,
                    )

                else:
                    trace = dict(type=trace_type, showlegend=False, **kwargs_trace)

                if x:
                    trace["x"] = group[x]
                if y:
                    trace["y"] = group[y]
                trace = _make_trace_for_scatter(
                    trace, trace_type, df[color_name], **kwargs_marker
                )

                fig.append_trace(trace, row_count + 1, col_count + 1)
                if row_count == 0:
                    label = _return_label(
                        col_values[col_count], facet_col_labels, facet_col
                    )
                    annotations.append(
                        _annotation_dict(
                            label,
                            col_count + 1,
                            num_of_cols,
                            SUBPLOT_SPACING,
                            row_col="col",
                            flipped=flipped_cols,
                        )
                    )
            label = _return_label(row_values[row_count], facet_row_labels, facet_row)
            annotations.append(
                _annotation_dict(
                    row_values[row_count],
                    num_of_rows - row_count,
                    num_of_rows,
                    SUBPLOT_SPACING,
                    row_col="row",
                    flipped=flipped_rows,
                )
            )

    return fig, annotations


def _facet_grid(
    df,
    x,
    y,
    facet_row,
    facet_col,
    num_of_rows,
    num_of_cols,
    facet_row_labels,
    facet_col_labels,
    trace_type,
    flipped_rows,
    flipped_cols,
    show_boxes,
    SUBPLOT_SPACING,
    marker_color,
    kwargs_trace,
    kwargs_marker,
):
    fig = make_subplots(
        rows=num_of_rows,
        cols=num_of_cols,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=SUBPLOT_SPACING,
        vertical_spacing=SUBPLOT_SPACING,
        print_grid=False,
    )
    annotations = []
    if not facet_row and not facet_col:
        trace = dict(
            type=trace_type,
            marker=dict(color=marker_color, line=kwargs_marker["line"]),
            **kwargs_trace,
        )

        if x:
            trace["x"] = df[x]
        if y:
            trace["y"] = df[y]
        trace = _make_trace_for_scatter(
            trace, trace_type, marker_color, **kwargs_marker
        )

        fig.append_trace(trace, 1, 1)

    elif (facet_row and not facet_col) or (not facet_row and facet_col):
        groups_by_facet = list(df.groupby(facet_row if facet_row else facet_col))
        for j, group in enumerate(groups_by_facet):
            trace = dict(
                type=trace_type,
                marker=dict(color=marker_color, line=kwargs_marker["line"]),
                **kwargs_trace,
            )

            if x:
                trace["x"] = group[1][x]
            if y:
                trace["y"] = group[1][y]
            trace = _make_trace_for_scatter(
                trace, trace_type, marker_color, **kwargs_marker
            )

            fig.append_trace(
                trace, j + 1 if facet_row else 1, 1 if facet_row else j + 1
            )

            label = _return_label(
                group[0],
                facet_row_labels if facet_row else facet_col_labels,
                facet_row if facet_row else facet_col,
            )

            annotations.append(
                _annotation_dict(
                    label,
                    num_of_rows - j if facet_row else j + 1,
                    num_of_rows if facet_row else num_of_cols,
                    SUBPLOT_SPACING,
                    "row" if facet_row else "col",
                    flipped_rows,
                )
            )

    elif facet_row and facet_col:
        groups_by_facets = list(df.groupby([facet_row, facet_col]))
        tuple_to_facet_group = {item[0]: item[1] for item in groups_by_facets}

        row_values = df[facet_row].unique()
        col_values = df[facet_col].unique()
        for row_count, x_val in enumerate(row_values):
            for col_count, y_val in enumerate(col_values):
                try:
                    group = tuple_to_facet_group[(x_val, y_val)]
                except KeyError:
                    group = pd.DataFrame([[None, None]], columns=[x, y])
                trace = dict(
                    type=trace_type,
                    marker=dict(color=marker_color, line=kwargs_marker["line"]),
                    **kwargs_trace,
                )
                if x:
                    trace["x"] = group[x]
                if y:
                    trace["y"] = group[y]
                trace = _make_trace_for_scatter(
                    trace, trace_type, marker_color, **kwargs_marker
                )

                fig.append_trace(trace, row_count + 1, col_count + 1)
                if row_count == 0:
                    label = _return_label(
                        col_values[col_count], facet_col_labels, facet_col
                    )
                    annotations.append(
                        _annotation_dict(
                            label,
                            col_count + 1,
                            num_of_cols,
                            SUBPLOT_SPACING,
                            row_col="col",
                            flipped=flipped_cols,
                        )
                    )

            label = _return_label(row_values[row_count], facet_row_labels, facet_row)
            annotations.append(
                _annotation_dict(
                    label,
                    num_of_rows - row_count,
                    num_of_rows,
                    SUBPLOT_SPACING,
                    row_col="row",
                    flipped=flipped_rows,
                )
            )

    return fig, annotations


def create_facet_grid(
    df,
    x=None,
    y=None,
    facet_row=None,
    facet_col=None,
    color_name=None,
    colormap=None,
    color_is_cat=False,
    facet_row_labels=None,
    facet_col_labels=None,
    height=None,
    width=None,
    trace_type="scatter",
    scales="fixed",
    dtick_x=None,
    dtick_y=None,
    show_boxes=True,
    ggplot2=False,
    binsize=1,
    **kwargs,
):
    """
    Returns figure for facet grid; **this function is deprecated**, since
    plotly.express functions should be used instead, for example

    >>> import plotly.express as px
    >>> tips = px.data.tips()
    >>> fig = px.scatter(tips,
    ...     x='total_bill',
    ...     y='tip',
    ...     facet_row='sex',
    ...     facet_col='smoker',
    ...     color='size')


    :param (pd.DataFrame) df: the dataframe of columns for the facet grid.
    :param (str) x: the name of the dataframe column for the x axis data.
    :param (str) y: the name of the dataframe column for the y axis data.
    :param (str) facet_row: the name of the dataframe column that is used to
        facet the grid into row panels.
    :param (str) facet_col: the name of the dataframe column that is used to
        facet the grid into column panels.
    :param (str) color_name: the name of your dataframe column that will
        function as the colormap variable.
    :param (str|list|dict) colormap: the param that determines how the
        color_name column colors the data. If the dataframe contains numeric
        data, then a dictionary of colors will group the data categorically
        while a Plotly Colorscale name or a custom colorscale will treat it
        numerically. To learn more about colors and types of colormap, run
        `help(plotly.colors)`.
    :param (bool) color_is_cat: determines whether a numerical column for the
        colormap will be treated as categorical (True) or sequential (False).
            Default = False.
    :param (str|dict) facet_row_labels: set to either 'name' or a dictionary
        of all the unique values in the faceting row mapped to some text to
        show up in the label annotations. If None, labeling works like usual.
    :param (str|dict) facet_col_labels: set to either 'name' or a dictionary
        of all the values in the faceting row mapped to some text to show up
        in the label annotations. If None, labeling works like usual.
    :param (int) height: the height of the facet grid figure.
    :param (int) width: the width of the facet grid figure.
    :param (str) trace_type: decides the type of plot to appear in the
        facet grid. The options are 'scatter', 'scattergl', 'histogram',
        'bar', and 'box'.
        Default = 'scatter'.
    :param (str) scales: determines if axes have fixed ranges or not. Valid
        settings are 'fixed' (all axes fixed), 'free_x' (x axis free only),
        'free_y' (y axis free only) or 'free' (both axes free).
    :param (float) dtick_x: determines the distance between each tick on the
        x-axis. Default is None which means dtick_x is set automatically.
    :param (float) dtick_y: determines the distance between each tick on the
        y-axis. Default is None which means dtick_y is set automatically.
    :param (bool) show_boxes: draws grey boxes behind the facet titles.
    :param (bool) ggplot2: draws the facet grid in the style of `ggplot2`. See
        http://ggplot2.tidyverse.org/reference/facet_grid.html for reference.
        Default = False
    :param (int) binsize: groups all data into bins of a given length.
    :param (dict) kwargs: a dictionary of scatterplot arguments.

    Examples 1: One Way Faceting

    >>> import plotly.figure_factory as ff
    >>> import pandas as pd
    >>> mpg = pd.read_table('https://raw.githubusercontent.com/plotly/datasets/master/mpg_2017.txt')

    >>> fig = ff.create_facet_grid(
    ...     mpg,
    ...     x='displ',
    ...     y='cty',
    ...     facet_col='cyl',
    ... )
    >>> fig.show()

    Example 2: Two Way Faceting

    >>> import plotly.figure_factory as ff

    >>> import pandas as pd

    >>> mpg = pd.read_table('https://raw.githubusercontent.com/plotly/datasets/master/mpg_2017.txt')

    >>> fig = ff.create_facet_grid(
    ...     mpg,
    ...     x='displ',
    ...     y='cty',
    ...     facet_row='drv',
    ...     facet_col='cyl',
    ... )
    >>> fig.show()

    Example 3: Categorical Coloring

    >>> import plotly.figure_factory as ff
    >>> import pandas as pd
    >>> mtcars = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/mtcars.csv')
    >>> mtcars.cyl = mtcars.cyl.astype(str)
    >>> fig = ff.create_facet_grid(
    ...     mtcars,
    ...     x='mpg',
    ...     y='wt',
    ...     facet_col='cyl',
    ...     color_name='cyl',
    ...     color_is_cat=True,
    ... )
    >>> fig.show()


    """
    if not pd:
        raise ImportError("'pandas' must be installed for this figure_factory.")

    if not isinstance(df, pd.DataFrame):
        raise exceptions.PlotlyError("You must input a pandas DataFrame.")

    # make sure all columns are of homogenous datatype
    utils.validate_dataframe(df)

    if trace_type in ["scatter", "scattergl"]:
        if not x or not y:
            raise exceptions.PlotlyError(
                "You need to input 'x' and 'y' if you are you are using a "
                "trace_type of 'scatter' or 'scattergl'."
            )

    for key in [x, y, facet_row, facet_col, color_name]:
        if key is not None:
            try:
                df[key]
            except KeyError:
                raise exceptions.PlotlyError(
                    "x, y, facet_row, facet_col and color_name must be keys "
                    "in your dataframe."
                )
    # autoscale histogram bars
    if trace_type not in ["scatter", "scattergl"]:
        scales = "free"

    # validate scales
    if scales not in ["fixed", "free_x", "free_y", "free"]:
        raise exceptions.PlotlyError(
            "'scales' must be set to 'fixed', 'free_x', 'free_y' and 'free'."
        )

    if trace_type not in VALID_TRACE_TYPES:
        raise exceptions.PlotlyError(
            "'trace_type' must be in {}".format(VALID_TRACE_TYPES)
        )

    if trace_type == "histogram":
        SUBPLOT_SPACING = 0.06
    else:
        SUBPLOT_SPACING = 0.015

    # seperate kwargs for marker and else
    if "marker" in kwargs:
        kwargs_marker = kwargs["marker"]
    else:
        kwargs_marker = {}
    marker_color = kwargs_marker.pop("color", None)
    kwargs.pop("marker", None)
    kwargs_trace = kwargs

    if "size" not in kwargs_marker:
        if ggplot2:
            kwargs_marker["size"] = 5
        else:
            kwargs_marker["size"] = 8

    if "opacity" not in kwargs_marker:
        if not ggplot2:
            kwargs_trace["opacity"] = 0.6

    if "line" not in kwargs_marker:
        if not ggplot2:
            kwargs_marker["line"] = {"color": "darkgrey", "width": 1}
        else:
            kwargs_marker["line"] = {}

    # default marker size
    if not ggplot2:
        if not marker_color:
            marker_color = "rgb(31, 119, 180)"
    else:
        marker_color = "rgb(0, 0, 0)"

    num_of_rows = 1
    num_of_cols = 1
    flipped_rows = False
    flipped_cols = False
    if facet_row:
        num_of_rows = len(df[facet_row].unique())
        flipped_rows = _is_flipped(num_of_rows)
        if isinstance(facet_row_labels, dict):
            for key in df[facet_row].unique():
                if key not in facet_row_labels.keys():
                    unique_keys = df[facet_row].unique().tolist()
                    raise exceptions.PlotlyError(CUSTOM_LABEL_ERROR.format(unique_keys))
    if facet_col:
        num_of_cols = len(df[facet_col].unique())
        flipped_cols = _is_flipped(num_of_cols)
        if isinstance(facet_col_labels, dict):
            for key in df[facet_col].unique():
                if key not in facet_col_labels.keys():
                    unique_keys = df[facet_col].unique().tolist()
                    raise exceptions.PlotlyError(CUSTOM_LABEL_ERROR.format(unique_keys))
    show_legend = False
    if color_name:
        if isinstance(df[color_name].iloc[0], str) or color_is_cat:
            show_legend = True
            if isinstance(colormap, dict):
                clrs.validate_colors_dict(colormap, "rgb")

                for val in df[color_name].unique():
                    if val not in colormap.keys():
                        raise exceptions.PlotlyError(
                            "If using 'colormap' as a dictionary, make sure "
                            "all the values of the colormap column are in "
                            "the keys of your dictionary."
                        )
            else:
                # use default plotly colors for dictionary
                default_colors = clrs.DEFAULT_PLOTLY_COLORS
                colormap = {}
                j = 0
                for val in df[color_name].unique():
                    if j >= len(default_colors):
                        j = 0
                    colormap[val] = default_colors[j]
                    j += 1
            fig, annotations = _facet_grid_color_categorical(
                df,
                x,
                y,
                facet_row,
                facet_col,
                color_name,
                colormap,
                num_of_rows,
                num_of_cols,
                facet_row_labels,
                facet_col_labels,
                trace_type,
                flipped_rows,
                flipped_cols,
                show_boxes,
                SUBPLOT_SPACING,
                marker_color,
                kwargs_trace,
                kwargs_marker,
            )

        elif isinstance(df[color_name].iloc[0], Number):
            if isinstance(colormap, dict):
                show_legend = True
                clrs.validate_colors_dict(colormap, "rgb")

                for val in df[color_name].unique():
                    if val not in colormap.keys():
                        raise exceptions.PlotlyError(
                            "If using 'colormap' as a dictionary, make sure "
                            "all the values of the colormap column are in "
                            "the keys of your dictionary."
                        )
                fig, annotations = _facet_grid_color_categorical(
                    df,
                    x,
                    y,
                    facet_row,
                    facet_col,
                    color_name,
                    colormap,
                    num_of_rows,
                    num_of_cols,
                    facet_row_labels,
                    facet_col_labels,
                    trace_type,
                    flipped_rows,
                    flipped_cols,
                    show_boxes,
                    SUBPLOT_SPACING,
                    marker_color,
                    kwargs_trace,
                    kwargs_marker,
                )

            elif isinstance(colormap, list):
                colorscale_list = colormap
                clrs.validate_colorscale(colorscale_list)

                fig, annotations = _facet_grid_color_numerical(
                    df,
                    x,
                    y,
                    facet_row,
                    facet_col,
                    color_name,
                    colorscale_list,
                    num_of_rows,
                    num_of_cols,
                    facet_row_labels,
                    facet_col_labels,
                    trace_type,
                    flipped_rows,
                    flipped_cols,
                    show_boxes,
                    SUBPLOT_SPACING,
                    marker_color,
                    kwargs_trace,
                    kwargs_marker,
                )
            elif isinstance(colormap, str):
                if colormap in clrs.PLOTLY_SCALES.keys():
                    colorscale_list = clrs.PLOTLY_SCALES[colormap]
                else:
                    raise exceptions.PlotlyError(
                        "If 'colormap' is a string, it must be the name "
                        "of a Plotly Colorscale. The available colorscale "
                        "names are {}".format(clrs.PLOTLY_SCALES.keys())
                    )
                fig, annotations = _facet_grid_color_numerical(
                    df,
                    x,
                    y,
                    facet_row,
                    facet_col,
                    color_name,
                    colorscale_list,
                    num_of_rows,
                    num_of_cols,
                    facet_row_labels,
                    facet_col_labels,
                    trace_type,
                    flipped_rows,
                    flipped_cols,
                    show_boxes,
                    SUBPLOT_SPACING,
                    marker_color,
                    kwargs_trace,
                    kwargs_marker,
                )
            else:
                colorscale_list = clrs.PLOTLY_SCALES["Reds"]
                fig, annotations = _facet_grid_color_numerical(
                    df,
                    x,
                    y,
                    facet_row,
                    facet_col,
                    color_name,
                    colorscale_list,
                    num_of_rows,
                    num_of_cols,
                    facet_row_labels,
                    facet_col_labels,
                    trace_type,
                    flipped_rows,
                    flipped_cols,
                    show_boxes,
                    SUBPLOT_SPACING,
                    marker_color,
                    kwargs_trace,
                    kwargs_marker,
                )

    else:
        fig, annotations = _facet_grid(
            df,
            x,
            y,
            facet_row,
            facet_col,
            num_of_rows,
            num_of_cols,
            facet_row_labels,
            facet_col_labels,
            trace_type,
            flipped_rows,
            flipped_cols,
            show_boxes,
            SUBPLOT_SPACING,
            marker_color,
            kwargs_trace,
            kwargs_marker,
        )

    if not height:
        height = max(600, 100 * num_of_rows)
    if not width:
        width = max(600, 100 * num_of_cols)

    fig["layout"].update(
        height=height, width=width, title="", paper_bgcolor="rgb(251, 251, 251)"
    )
    if ggplot2:
        fig["layout"].update(
            plot_bgcolor=PLOT_BGCOLOR,
            paper_bgcolor="rgb(255, 255, 255)",
            hovermode="closest",
        )

    # axis titles
    x_title_annot = _axis_title_annotation(x, "x")
    y_title_annot = _axis_title_annotation(y, "y")

    # annotations
    annotations.append(x_title_annot)
    annotations.append(y_title_annot)

    # legend
    fig["layout"]["showlegend"] = show_legend
    fig["layout"]["legend"]["bgcolor"] = LEGEND_COLOR
    fig["layout"]["legend"]["borderwidth"] = LEGEND_BORDER_WIDTH
    fig["layout"]["legend"]["x"] = 1.05
    fig["layout"]["legend"]["y"] = 1
    fig["layout"]["legend"]["yanchor"] = "top"

    if show_legend:
        fig["layout"]["showlegend"] = show_legend
        if ggplot2:
            if color_name:
                legend_annot = _legend_annotation(color_name)
                annotations.append(legend_annot)
            fig["layout"]["margin"]["r"] = 150

    # assign annotations to figure
    fig["layout"]["annotations"] = annotations

    # add shaded boxes behind axis titles
    if show_boxes and ggplot2:
        _add_shapes_to_fig(fig, ANNOT_RECT_COLOR, flipped_rows, flipped_cols)

    # all xaxis and yaxis labels
    axis_labels = {"x": [], "y": []}
    for key in fig["layout"]:
        if "xaxis" in key:
            axis_labels["x"].append(key)
        elif "yaxis" in key:
            axis_labels["y"].append(key)

    string_number_in_data = False
    for var in [v for v in [x, y] if v]:
        if isinstance(df[var].tolist()[0], str):
            for item in df[var]:
                try:
                    int(item)
                    string_number_in_data = True
                except ValueError:
                    pass

    if string_number_in_data:
        for x_y in axis_labels.keys():
            for axis_name in axis_labels[x_y]:
                fig["layout"][axis_name]["type"] = "category"

    if scales == "fixed":
        fixed_axes = ["x", "y"]
    elif scales == "free_x":
        fixed_axes = ["y"]
    elif scales == "free_y":
        fixed_axes = ["x"]
    elif scales == "free":
        fixed_axes = []

    # fixed ranges
    for x_y in fixed_axes:
        min_ranges = []
        max_ranges = []
        for trace in fig["data"]:
            if trace[x_y] is not None and len(trace[x_y]) > 0:
                min_ranges.append(min(trace[x_y]))
                max_ranges.append(max(trace[x_y]))
        while None in min_ranges:
            min_ranges.remove(None)
        while None in max_ranges:
            max_ranges.remove(None)

        min_range = min(min_ranges)
        max_range = max(max_ranges)

        range_are_numbers = isinstance(min_range, Number) and isinstance(
            max_range, Number
        )

        if range_are_numbers:
            min_range = math.floor(min_range)
            max_range = math.ceil(max_range)

            # extend widen frame by 5% on each side
            min_range -= 0.05 * (max_range - min_range)
            max_range += 0.05 * (max_range - min_range)

            if x_y == "x":
                if dtick_x:
                    dtick = dtick_x
                else:
                    dtick = math.floor((max_range - min_range) / MAX_TICKS_PER_AXIS)
            elif x_y == "y":
                if dtick_y:
                    dtick = dtick_y
                else:
                    dtick = math.floor((max_range - min_range) / MAX_TICKS_PER_AXIS)
        else:
            dtick = 1

        for axis_title in axis_labels[x_y]:
            fig["layout"][axis_title]["dtick"] = dtick
            fig["layout"][axis_title]["ticklen"] = 0
            fig["layout"][axis_title]["zeroline"] = False
            if ggplot2:
                fig["layout"][axis_title]["tickwidth"] = 1
                fig["layout"][axis_title]["ticklen"] = 4
                fig["layout"][axis_title]["gridwidth"] = GRID_WIDTH

                fig["layout"][axis_title]["gridcolor"] = GRID_COLOR
                fig["layout"][axis_title]["gridwidth"] = 2
                fig["layout"][axis_title]["tickfont"] = {
                    "color": TICK_COLOR,
                    "size": 10,
                }

        # insert ranges into fig
        if x_y in fixed_axes:
            for key in fig["layout"]:
                if "{}axis".format(x_y) in key and range_are_numbers:
                    fig["layout"][key]["range"] = [min_range, max_range]

    return fig
