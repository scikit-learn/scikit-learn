from plotly import exceptions, optional_imports
import plotly.colors as clrs
from plotly.figure_factory import utils
from plotly.graph_objs import graph_objs
from plotly.subplots import make_subplots

pd = optional_imports.get_module("pandas")

DIAG_CHOICES = ["scatter", "histogram", "box"]
VALID_COLORMAP_TYPES = ["cat", "seq"]


def endpts_to_intervals(endpts):
    """
    Returns a list of intervals for categorical colormaps

    Accepts a list or tuple of sequentially increasing numbers and returns
    a list representation of the mathematical intervals with these numbers
    as endpoints. For example, [1, 6] returns [[-inf, 1], [1, 6], [6, inf]]

    :raises: (PlotlyError) If input is not a list or tuple
    :raises: (PlotlyError) If the input contains a string
    :raises: (PlotlyError) If any number does not increase after the
        previous one in the sequence
    """
    length = len(endpts)
    # Check if endpts is a list or tuple
    if not (isinstance(endpts, (tuple)) or isinstance(endpts, (list))):
        raise exceptions.PlotlyError(
            "The intervals_endpts argument must "
            "be a list or tuple of a sequence "
            "of increasing numbers."
        )
    # Check if endpts contains only numbers
    for item in endpts:
        if isinstance(item, str):
            raise exceptions.PlotlyError(
                "The intervals_endpts argument "
                "must be a list or tuple of a "
                "sequence of increasing "
                "numbers."
            )
    # Check if numbers in endpts are increasing
    for k in range(length - 1):
        if endpts[k] >= endpts[k + 1]:
            raise exceptions.PlotlyError(
                "The intervals_endpts argument "
                "must be a list or tuple of a "
                "sequence of increasing "
                "numbers."
            )
    else:
        intervals = []
        # add -inf to intervals
        intervals.append([float("-inf"), endpts[0]])
        for k in range(length - 1):
            interval = []
            interval.append(endpts[k])
            interval.append(endpts[k + 1])
            intervals.append(interval)
        # add +inf to intervals
        intervals.append([endpts[length - 1], float("inf")])
        return intervals


def hide_tick_labels_from_box_subplots(fig):
    """
    Hides tick labels for box plots in scatterplotmatrix subplots.
    """
    boxplot_xaxes = []
    for trace in fig["data"]:
        if trace["type"] == "box":
            # stores the xaxes which correspond to boxplot subplots
            # since we use xaxis1, xaxis2, etc, in plotly.py
            boxplot_xaxes.append("xaxis{}".format(trace["xaxis"][1:]))
    for xaxis in boxplot_xaxes:
        fig["layout"][xaxis]["showticklabels"] = False


def validate_scatterplotmatrix(df, index, diag, colormap_type, **kwargs):
    """
    Validates basic inputs for FigureFactory.create_scatterplotmatrix()

    :raises: (PlotlyError) If pandas is not imported
    :raises: (PlotlyError) If pandas dataframe is not inputted
    :raises: (PlotlyError) If pandas dataframe has <= 1 columns
    :raises: (PlotlyError) If diagonal plot choice (diag) is not one of
        the viable options
    :raises: (PlotlyError) If colormap_type is not a valid choice
    :raises: (PlotlyError) If kwargs contains 'size', 'color' or
        'colorscale'
    """
    if not pd:
        raise ImportError(
            "FigureFactory.scatterplotmatrix requires a pandas DataFrame."
        )

    # Check if pandas dataframe
    if not isinstance(df, pd.core.frame.DataFrame):
        raise exceptions.PlotlyError(
            "Dataframe not inputed. Please "
            "use a pandas dataframe to pro"
            "duce a scatterplot matrix."
        )

    # Check if dataframe is 1 column or less
    if len(df.columns) <= 1:
        raise exceptions.PlotlyError(
            "Dataframe has only one column. To "
            "use the scatterplot matrix, use at "
            "least 2 columns."
        )

    # Check that diag parameter is a valid selection
    if diag not in DIAG_CHOICES:
        raise exceptions.PlotlyError(
            "Make sure diag is set to one of {}".format(DIAG_CHOICES)
        )

    # Check that colormap_types is a valid selection
    if colormap_type not in VALID_COLORMAP_TYPES:
        raise exceptions.PlotlyError(
            "Must choose a valid colormap type. "
            "Either 'cat' or 'seq' for a cate"
            "gorical and sequential colormap "
            "respectively."
        )

    # Check for not 'size' or 'color' in 'marker' of **kwargs
    if "marker" in kwargs:
        FORBIDDEN_PARAMS = ["size", "color", "colorscale"]
        if any(param in kwargs["marker"] for param in FORBIDDEN_PARAMS):
            raise exceptions.PlotlyError(
                "Your kwargs dictionary cannot "
                "include the 'size', 'color' or "
                "'colorscale' key words inside "
                "the marker dict since 'size' is "
                "already an argument of the "
                "scatterplot matrix function and "
                "both 'color' and 'colorscale "
                "are set internally."
            )


def scatterplot(dataframe, headers, diag, size, height, width, title, **kwargs):
    """
    Refer to FigureFactory.create_scatterplotmatrix() for docstring

    Returns fig for scatterplotmatrix without index

    """
    dim = len(dataframe)
    fig = make_subplots(rows=dim, cols=dim, print_grid=False)
    trace_list = []
    # Insert traces into trace_list
    for listy in dataframe:
        for listx in dataframe:
            if (listx == listy) and (diag == "histogram"):
                trace = graph_objs.Histogram(x=listx, showlegend=False)
            elif (listx == listy) and (diag == "box"):
                trace = graph_objs.Box(y=listx, name=None, showlegend=False)
            else:
                if "marker" in kwargs:
                    kwargs["marker"]["size"] = size
                    trace = graph_objs.Scatter(
                        x=listx, y=listy, mode="markers", showlegend=False, **kwargs
                    )
                    trace_list.append(trace)
                else:
                    trace = graph_objs.Scatter(
                        x=listx,
                        y=listy,
                        mode="markers",
                        marker=dict(size=size),
                        showlegend=False,
                        **kwargs,
                    )
            trace_list.append(trace)

    trace_index = 0
    indices = range(1, dim + 1)
    for y_index in indices:
        for x_index in indices:
            fig.append_trace(trace_list[trace_index], y_index, x_index)
            trace_index += 1

    # Insert headers into the figure
    for j in range(dim):
        xaxis_key = "xaxis{}".format((dim * dim) - dim + 1 + j)
        fig["layout"][xaxis_key].update(title=headers[j])
    for j in range(dim):
        yaxis_key = "yaxis{}".format(1 + (dim * j))
        fig["layout"][yaxis_key].update(title=headers[j])

    fig["layout"].update(height=height, width=width, title=title, showlegend=True)

    hide_tick_labels_from_box_subplots(fig)

    return fig


def scatterplot_dict(
    dataframe,
    headers,
    diag,
    size,
    height,
    width,
    title,
    index,
    index_vals,
    endpts,
    colormap,
    colormap_type,
    **kwargs,
):
    """
    Refer to FigureFactory.create_scatterplotmatrix() for docstring

    Returns fig for scatterplotmatrix with both index and colormap picked.
    Used if colormap is a dictionary with index values as keys pointing to
    colors. Forces colormap_type to behave categorically because it would
    not make sense colors are assigned to each index value and thus
    implies that a categorical approach should be taken

    """

    theme = colormap
    dim = len(dataframe)
    fig = make_subplots(rows=dim, cols=dim, print_grid=False)
    trace_list = []
    legend_param = 0
    # Work over all permutations of list pairs
    for listy in dataframe:
        for listx in dataframe:
            # create a dictionary for index_vals
            unique_index_vals = {}
            for name in index_vals:
                if name not in unique_index_vals:
                    unique_index_vals[name] = []

            # Fill all the rest of the names into the dictionary
            for name in sorted(unique_index_vals.keys()):
                new_listx = []
                new_listy = []
                for j in range(len(index_vals)):
                    if index_vals[j] == name:
                        new_listx.append(listx[j])
                        new_listy.append(listy[j])
                # Generate trace with VISIBLE icon
                if legend_param == 1:
                    if (listx == listy) and (diag == "histogram"):
                        trace = graph_objs.Histogram(
                            x=new_listx, marker=dict(color=theme[name]), showlegend=True
                        )
                    elif (listx == listy) and (diag == "box"):
                        trace = graph_objs.Box(
                            y=new_listx,
                            name=None,
                            marker=dict(color=theme[name]),
                            showlegend=True,
                        )
                    else:
                        if "marker" in kwargs:
                            kwargs["marker"]["size"] = size
                            kwargs["marker"]["color"] = theme[name]
                            trace = graph_objs.Scatter(
                                x=new_listx,
                                y=new_listy,
                                mode="markers",
                                name=name,
                                showlegend=True,
                                **kwargs,
                            )
                        else:
                            trace = graph_objs.Scatter(
                                x=new_listx,
                                y=new_listy,
                                mode="markers",
                                name=name,
                                marker=dict(size=size, color=theme[name]),
                                showlegend=True,
                                **kwargs,
                            )
                # Generate trace with INVISIBLE icon
                else:
                    if (listx == listy) and (diag == "histogram"):
                        trace = graph_objs.Histogram(
                            x=new_listx,
                            marker=dict(color=theme[name]),
                            showlegend=False,
                        )
                    elif (listx == listy) and (diag == "box"):
                        trace = graph_objs.Box(
                            y=new_listx,
                            name=None,
                            marker=dict(color=theme[name]),
                            showlegend=False,
                        )
                    else:
                        if "marker" in kwargs:
                            kwargs["marker"]["size"] = size
                            kwargs["marker"]["color"] = theme[name]
                            trace = graph_objs.Scatter(
                                x=new_listx,
                                y=new_listy,
                                mode="markers",
                                name=name,
                                showlegend=False,
                                **kwargs,
                            )
                        else:
                            trace = graph_objs.Scatter(
                                x=new_listx,
                                y=new_listy,
                                mode="markers",
                                name=name,
                                marker=dict(size=size, color=theme[name]),
                                showlegend=False,
                                **kwargs,
                            )
                # Push the trace into dictionary
                unique_index_vals[name] = trace
            trace_list.append(unique_index_vals)
            legend_param += 1

    trace_index = 0
    indices = range(1, dim + 1)
    for y_index in indices:
        for x_index in indices:
            for name in sorted(trace_list[trace_index].keys()):
                fig.append_trace(trace_list[trace_index][name], y_index, x_index)
            trace_index += 1

    # Insert headers into the figure
    for j in range(dim):
        xaxis_key = "xaxis{}".format((dim * dim) - dim + 1 + j)
        fig["layout"][xaxis_key].update(title=headers[j])

    for j in range(dim):
        yaxis_key = "yaxis{}".format(1 + (dim * j))
        fig["layout"][yaxis_key].update(title=headers[j])

    hide_tick_labels_from_box_subplots(fig)

    if diag == "histogram":
        fig["layout"].update(
            height=height, width=width, title=title, showlegend=True, barmode="stack"
        )
        return fig

    else:
        fig["layout"].update(height=height, width=width, title=title, showlegend=True)
        return fig


def scatterplot_theme(
    dataframe,
    headers,
    diag,
    size,
    height,
    width,
    title,
    index,
    index_vals,
    endpts,
    colormap,
    colormap_type,
    **kwargs,
):
    """
    Refer to FigureFactory.create_scatterplotmatrix() for docstring

    Returns fig for scatterplotmatrix with both index and colormap picked

    """

    # Check if index is made of string values
    if isinstance(index_vals[0], str):
        unique_index_vals = []
        for name in index_vals:
            if name not in unique_index_vals:
                unique_index_vals.append(name)
        n_colors_len = len(unique_index_vals)

        # Convert colormap to list of n RGB tuples
        if colormap_type == "seq":
            foo = clrs.color_parser(colormap, clrs.unlabel_rgb)
            foo = clrs.n_colors(foo[0], foo[1], n_colors_len)
            theme = clrs.color_parser(foo, clrs.label_rgb)

        if colormap_type == "cat":
            # leave list of colors the same way
            theme = colormap

        dim = len(dataframe)
        fig = make_subplots(rows=dim, cols=dim, print_grid=False)
        trace_list = []
        legend_param = 0
        # Work over all permutations of list pairs
        for listy in dataframe:
            for listx in dataframe:
                # create a dictionary for index_vals
                unique_index_vals = {}
                for name in index_vals:
                    if name not in unique_index_vals:
                        unique_index_vals[name] = []

                c_indx = 0  # color index
                # Fill all the rest of the names into the dictionary
                for name in sorted(unique_index_vals.keys()):
                    new_listx = []
                    new_listy = []
                    for j in range(len(index_vals)):
                        if index_vals[j] == name:
                            new_listx.append(listx[j])
                            new_listy.append(listy[j])
                    # Generate trace with VISIBLE icon
                    if legend_param == 1:
                        if (listx == listy) and (diag == "histogram"):
                            trace = graph_objs.Histogram(
                                x=new_listx,
                                marker=dict(color=theme[c_indx]),
                                showlegend=True,
                            )
                        elif (listx == listy) and (diag == "box"):
                            trace = graph_objs.Box(
                                y=new_listx,
                                name=None,
                                marker=dict(color=theme[c_indx]),
                                showlegend=True,
                            )
                        else:
                            if "marker" in kwargs:
                                kwargs["marker"]["size"] = size
                                kwargs["marker"]["color"] = theme[c_indx]
                                trace = graph_objs.Scatter(
                                    x=new_listx,
                                    y=new_listy,
                                    mode="markers",
                                    name=name,
                                    showlegend=True,
                                    **kwargs,
                                )
                            else:
                                trace = graph_objs.Scatter(
                                    x=new_listx,
                                    y=new_listy,
                                    mode="markers",
                                    name=name,
                                    marker=dict(size=size, color=theme[c_indx]),
                                    showlegend=True,
                                    **kwargs,
                                )
                    # Generate trace with INVISIBLE icon
                    else:
                        if (listx == listy) and (diag == "histogram"):
                            trace = graph_objs.Histogram(
                                x=new_listx,
                                marker=dict(color=theme[c_indx]),
                                showlegend=False,
                            )
                        elif (listx == listy) and (diag == "box"):
                            trace = graph_objs.Box(
                                y=new_listx,
                                name=None,
                                marker=dict(color=theme[c_indx]),
                                showlegend=False,
                            )
                        else:
                            if "marker" in kwargs:
                                kwargs["marker"]["size"] = size
                                kwargs["marker"]["color"] = theme[c_indx]
                                trace = graph_objs.Scatter(
                                    x=new_listx,
                                    y=new_listy,
                                    mode="markers",
                                    name=name,
                                    showlegend=False,
                                    **kwargs,
                                )
                            else:
                                trace = graph_objs.Scatter(
                                    x=new_listx,
                                    y=new_listy,
                                    mode="markers",
                                    name=name,
                                    marker=dict(size=size, color=theme[c_indx]),
                                    showlegend=False,
                                    **kwargs,
                                )
                    # Push the trace into dictionary
                    unique_index_vals[name] = trace
                    if c_indx >= (len(theme) - 1):
                        c_indx = -1
                    c_indx += 1
                trace_list.append(unique_index_vals)
                legend_param += 1

        trace_index = 0
        indices = range(1, dim + 1)
        for y_index in indices:
            for x_index in indices:
                for name in sorted(trace_list[trace_index].keys()):
                    fig.append_trace(trace_list[trace_index][name], y_index, x_index)
                trace_index += 1

        # Insert headers into the figure
        for j in range(dim):
            xaxis_key = "xaxis{}".format((dim * dim) - dim + 1 + j)
            fig["layout"][xaxis_key].update(title=headers[j])

        for j in range(dim):
            yaxis_key = "yaxis{}".format(1 + (dim * j))
            fig["layout"][yaxis_key].update(title=headers[j])

        hide_tick_labels_from_box_subplots(fig)

        if diag == "histogram":
            fig["layout"].update(
                height=height,
                width=width,
                title=title,
                showlegend=True,
                barmode="stack",
            )
            return fig

        elif diag == "box":
            fig["layout"].update(
                height=height, width=width, title=title, showlegend=True
            )
            return fig

        else:
            fig["layout"].update(
                height=height, width=width, title=title, showlegend=True
            )
            return fig

    else:
        if endpts:
            intervals = utils.endpts_to_intervals(endpts)

            # Convert colormap to list of n RGB tuples
            if colormap_type == "seq":
                foo = clrs.color_parser(colormap, clrs.unlabel_rgb)
                foo = clrs.n_colors(foo[0], foo[1], len(intervals))
                theme = clrs.color_parser(foo, clrs.label_rgb)

            if colormap_type == "cat":
                # leave list of colors the same way
                theme = colormap

            dim = len(dataframe)
            fig = make_subplots(rows=dim, cols=dim, print_grid=False)
            trace_list = []
            legend_param = 0
            # Work over all permutations of list pairs
            for listy in dataframe:
                for listx in dataframe:
                    interval_labels = {}
                    for interval in intervals:
                        interval_labels[str(interval)] = []

                    c_indx = 0  # color index
                    # Fill all the rest of the names into the dictionary
                    for interval in intervals:
                        new_listx = []
                        new_listy = []
                        for j in range(len(index_vals)):
                            if interval[0] < index_vals[j] <= interval[1]:
                                new_listx.append(listx[j])
                                new_listy.append(listy[j])
                        # Generate trace with VISIBLE icon
                        if legend_param == 1:
                            if (listx == listy) and (diag == "histogram"):
                                trace = graph_objs.Histogram(
                                    x=new_listx,
                                    marker=dict(color=theme[c_indx]),
                                    showlegend=True,
                                )
                            elif (listx == listy) and (diag == "box"):
                                trace = graph_objs.Box(
                                    y=new_listx,
                                    name=None,
                                    marker=dict(color=theme[c_indx]),
                                    showlegend=True,
                                )
                            else:
                                if "marker" in kwargs:
                                    kwargs["marker"]["size"] = size
                                    (kwargs["marker"]["color"]) = theme[c_indx]
                                    trace = graph_objs.Scatter(
                                        x=new_listx,
                                        y=new_listy,
                                        mode="markers",
                                        name=str(interval),
                                        showlegend=True,
                                        **kwargs,
                                    )
                                else:
                                    trace = graph_objs.Scatter(
                                        x=new_listx,
                                        y=new_listy,
                                        mode="markers",
                                        name=str(interval),
                                        marker=dict(size=size, color=theme[c_indx]),
                                        showlegend=True,
                                        **kwargs,
                                    )
                        # Generate trace with INVISIBLE icon
                        else:
                            if (listx == listy) and (diag == "histogram"):
                                trace = graph_objs.Histogram(
                                    x=new_listx,
                                    marker=dict(color=theme[c_indx]),
                                    showlegend=False,
                                )
                            elif (listx == listy) and (diag == "box"):
                                trace = graph_objs.Box(
                                    y=new_listx,
                                    name=None,
                                    marker=dict(color=theme[c_indx]),
                                    showlegend=False,
                                )
                            else:
                                if "marker" in kwargs:
                                    kwargs["marker"]["size"] = size
                                    (kwargs["marker"]["color"]) = theme[c_indx]
                                    trace = graph_objs.Scatter(
                                        x=new_listx,
                                        y=new_listy,
                                        mode="markers",
                                        name=str(interval),
                                        showlegend=False,
                                        **kwargs,
                                    )
                                else:
                                    trace = graph_objs.Scatter(
                                        x=new_listx,
                                        y=new_listy,
                                        mode="markers",
                                        name=str(interval),
                                        marker=dict(size=size, color=theme[c_indx]),
                                        showlegend=False,
                                        **kwargs,
                                    )
                        # Push the trace into dictionary
                        interval_labels[str(interval)] = trace
                        if c_indx >= (len(theme) - 1):
                            c_indx = -1
                        c_indx += 1
                    trace_list.append(interval_labels)
                    legend_param += 1

            trace_index = 0
            indices = range(1, dim + 1)
            for y_index in indices:
                for x_index in indices:
                    for interval in intervals:
                        fig.append_trace(
                            trace_list[trace_index][str(interval)], y_index, x_index
                        )
                    trace_index += 1

            # Insert headers into the figure
            for j in range(dim):
                xaxis_key = "xaxis{}".format((dim * dim) - dim + 1 + j)
                fig["layout"][xaxis_key].update(title=headers[j])
            for j in range(dim):
                yaxis_key = "yaxis{}".format(1 + (dim * j))
                fig["layout"][yaxis_key].update(title=headers[j])

            hide_tick_labels_from_box_subplots(fig)

            if diag == "histogram":
                fig["layout"].update(
                    height=height,
                    width=width,
                    title=title,
                    showlegend=True,
                    barmode="stack",
                )
                return fig

            elif diag == "box":
                fig["layout"].update(
                    height=height, width=width, title=title, showlegend=True
                )
                return fig

            else:
                fig["layout"].update(
                    height=height, width=width, title=title, showlegend=True
                )
                return fig

        else:
            theme = colormap

            # add a copy of rgb color to theme if it contains one color
            if len(theme) <= 1:
                theme.append(theme[0])

            color = []
            for incr in range(len(theme)):
                color.append([1.0 / (len(theme) - 1) * incr, theme[incr]])

            dim = len(dataframe)
            fig = make_subplots(rows=dim, cols=dim, print_grid=False)
            trace_list = []
            legend_param = 0
            # Run through all permutations of list pairs
            for listy in dataframe:
                for listx in dataframe:
                    # Generate trace with VISIBLE icon
                    if legend_param == 1:
                        if (listx == listy) and (diag == "histogram"):
                            trace = graph_objs.Histogram(
                                x=listx, marker=dict(color=theme[0]), showlegend=False
                            )
                        elif (listx == listy) and (diag == "box"):
                            trace = graph_objs.Box(
                                y=listx, marker=dict(color=theme[0]), showlegend=False
                            )
                        else:
                            if "marker" in kwargs:
                                kwargs["marker"]["size"] = size
                                kwargs["marker"]["color"] = index_vals
                                kwargs["marker"]["colorscale"] = color
                                kwargs["marker"]["showscale"] = True
                                trace = graph_objs.Scatter(
                                    x=listx,
                                    y=listy,
                                    mode="markers",
                                    showlegend=False,
                                    **kwargs,
                                )
                            else:
                                trace = graph_objs.Scatter(
                                    x=listx,
                                    y=listy,
                                    mode="markers",
                                    marker=dict(
                                        size=size,
                                        color=index_vals,
                                        colorscale=color,
                                        showscale=True,
                                    ),
                                    showlegend=False,
                                    **kwargs,
                                )
                    # Generate trace with INVISIBLE icon
                    else:
                        if (listx == listy) and (diag == "histogram"):
                            trace = graph_objs.Histogram(
                                x=listx, marker=dict(color=theme[0]), showlegend=False
                            )
                        elif (listx == listy) and (diag == "box"):
                            trace = graph_objs.Box(
                                y=listx, marker=dict(color=theme[0]), showlegend=False
                            )
                        else:
                            if "marker" in kwargs:
                                kwargs["marker"]["size"] = size
                                kwargs["marker"]["color"] = index_vals
                                kwargs["marker"]["colorscale"] = color
                                kwargs["marker"]["showscale"] = False
                                trace = graph_objs.Scatter(
                                    x=listx,
                                    y=listy,
                                    mode="markers",
                                    showlegend=False,
                                    **kwargs,
                                )
                            else:
                                trace = graph_objs.Scatter(
                                    x=listx,
                                    y=listy,
                                    mode="markers",
                                    marker=dict(
                                        size=size,
                                        color=index_vals,
                                        colorscale=color,
                                        showscale=False,
                                    ),
                                    showlegend=False,
                                    **kwargs,
                                )
                    # Push the trace into list
                    trace_list.append(trace)
                    legend_param += 1

            trace_index = 0
            indices = range(1, dim + 1)
            for y_index in indices:
                for x_index in indices:
                    fig.append_trace(trace_list[trace_index], y_index, x_index)
                    trace_index += 1

            # Insert headers into the figure
            for j in range(dim):
                xaxis_key = "xaxis{}".format((dim * dim) - dim + 1 + j)
                fig["layout"][xaxis_key].update(title=headers[j])
            for j in range(dim):
                yaxis_key = "yaxis{}".format(1 + (dim * j))
                fig["layout"][yaxis_key].update(title=headers[j])

            hide_tick_labels_from_box_subplots(fig)

            if diag == "histogram":
                fig["layout"].update(
                    height=height,
                    width=width,
                    title=title,
                    showlegend=True,
                    barmode="stack",
                )
                return fig

            elif diag == "box":
                fig["layout"].update(
                    height=height, width=width, title=title, showlegend=True
                )
                return fig

            else:
                fig["layout"].update(
                    height=height, width=width, title=title, showlegend=True
                )
                return fig


def create_scatterplotmatrix(
    df,
    index=None,
    endpts=None,
    diag="scatter",
    height=500,
    width=500,
    size=6,
    title="Scatterplot Matrix",
    colormap=None,
    colormap_type="cat",
    dataframe=None,
    headers=None,
    index_vals=None,
    **kwargs,
):
    """
    Returns data for a scatterplot matrix;
    **deprecated**,
    use instead the plotly.graph_objects trace
    :class:`plotly.graph_objects.Splom`.

    :param (array) df: array of the data with column headers
    :param (str) index: name of the index column in data array
    :param (list|tuple) endpts: takes an increasing sequece of numbers
        that defines intervals on the real line. They are used to group
        the entries in an index of numbers into their corresponding
        interval and therefore can be treated as categorical data
    :param (str) diag: sets the chart type for the main diagonal plots.
        The options are 'scatter', 'histogram' and 'box'.
    :param (int|float) height: sets the height of the chart
    :param (int|float) width: sets the width of the chart
    :param (float) size: sets the marker size (in px)
    :param (str) title: the title label of the scatterplot matrix
    :param (str|tuple|list|dict) colormap: either a plotly scale name,
        an rgb or hex color, a color tuple, a list of colors or a
        dictionary. An rgb color is of the form 'rgb(x, y, z)' where
        x, y and z belong to the interval [0, 255] and a color tuple is a
        tuple of the form (a, b, c) where a, b and c belong to [0, 1].
        If colormap is a list, it must contain valid color types as its
        members.
        If colormap is a dictionary, all the string entries in
        the index column must be a key in colormap. In this case, the
        colormap_type is forced to 'cat' or categorical
    :param (str) colormap_type: determines how colormap is interpreted.
        Valid choices are 'seq' (sequential) and 'cat' (categorical). If
        'seq' is selected, only the first two colors in colormap will be
        considered (when colormap is a list) and the index values will be
        linearly interpolated between those two colors. This option is
        forced if all index values are numeric.
        If 'cat' is selected, a color from colormap will be assigned to
        each category from index, including the intervals if endpts is
        being used
    :param (dict) **kwargs: a dictionary of scatterplot arguments
        The only forbidden parameters are 'size', 'color' and
        'colorscale' in 'marker'

    Example 1: Vanilla Scatterplot Matrix

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd

    >>> # Create dataframe
    >>> df = pd.DataFrame(np.random.randn(10, 2),
    ...                 columns=['Column 1', 'Column 2'])

    >>> # Create scatterplot matrix
    >>> fig = create_scatterplotmatrix(df)
    >>> fig.show()


    Example 2: Indexing a Column

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd

    >>> # Create dataframe with index
    >>> df = pd.DataFrame(np.random.randn(10, 2),
    ...                    columns=['A', 'B'])

    >>> # Add another column of strings to the dataframe
    >>> df['Fruit'] = pd.Series(['apple', 'apple', 'grape', 'apple', 'apple',
    ...                          'grape', 'pear', 'pear', 'apple', 'pear'])

    >>> # Create scatterplot matrix
    >>> fig = create_scatterplotmatrix(df, index='Fruit', size=10)
    >>> fig.show()


    Example 3: Styling the Diagonal Subplots

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd

    >>> # Create dataframe with index
    >>> df = pd.DataFrame(np.random.randn(10, 4),
    ...                    columns=['A', 'B', 'C', 'D'])

    >>> # Add another column of strings to the dataframe
    >>> df['Fruit'] = pd.Series(['apple', 'apple', 'grape', 'apple', 'apple',
    ...                          'grape', 'pear', 'pear', 'apple', 'pear'])

    >>> # Create scatterplot matrix
    >>> fig = create_scatterplotmatrix(df, diag='box', index='Fruit', height=1000,
    ...                                width=1000)
    >>> fig.show()


    Example 4: Use a Theme to Style the Subplots

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd

    >>> # Create dataframe with random data
    >>> df = pd.DataFrame(np.random.randn(100, 3),
    ...                    columns=['A', 'B', 'C'])

    >>> # Create scatterplot matrix using a built-in
    >>> # Plotly palette scale and indexing column 'A'
    >>> fig = create_scatterplotmatrix(df, diag='histogram', index='A',
    ...                                colormap='Blues', height=800, width=800)
    >>> fig.show()


    Example 5: Example 4 with Interval Factoring

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd

    >>> # Create dataframe with random data
    >>> df = pd.DataFrame(np.random.randn(100, 3),
    ...                    columns=['A', 'B', 'C'])

    >>> # Create scatterplot matrix using a list of 2 rgb tuples
    >>> # and endpoints at -1, 0 and 1
    >>> fig = create_scatterplotmatrix(df, diag='histogram', index='A',
    ...                                colormap=['rgb(140, 255, 50)',
    ...                                          'rgb(170, 60, 115)', '#6c4774',
    ...                                          (0.5, 0.1, 0.8)],
    ...                                endpts=[-1, 0, 1], height=800, width=800)
    >>> fig.show()


    Example 6: Using the colormap as a Dictionary

    >>> from plotly.graph_objs import graph_objs
    >>> from plotly.figure_factory import create_scatterplotmatrix

    >>> import numpy as np
    >>> import pandas as pd
    >>> import random

    >>> # Create dataframe with random data
    >>> df = pd.DataFrame(np.random.randn(100, 3),
    ...                    columns=['Column A',
    ...                             'Column B',
    ...                             'Column C'])

    >>> # Add new color column to dataframe
    >>> new_column = []
    >>> strange_colors = ['turquoise', 'limegreen', 'goldenrod']

    >>> for j in range(100):
    ...     new_column.append(random.choice(strange_colors))
    >>> df['Colors'] = pd.Series(new_column, index=df.index)

    >>> # Create scatterplot matrix using a dictionary of hex color values
    >>> # which correspond to actual color names in 'Colors' column
    >>> fig = create_scatterplotmatrix(
    ...     df, diag='box', index='Colors',
    ...     colormap= dict(
    ...         turquoise = '#00F5FF',
    ...         limegreen = '#32CD32',
    ...         goldenrod = '#DAA520'
    ...     ),
    ...     colormap_type='cat',
    ...     height=800, width=800
    ... )
    >>> fig.show()
    """
    # TODO: protected until #282
    if dataframe is None:
        dataframe = []
    if headers is None:
        headers = []
    if index_vals is None:
        index_vals = []

    validate_scatterplotmatrix(df, index, diag, colormap_type, **kwargs)

    # Validate colormap
    if isinstance(colormap, dict):
        colormap = clrs.validate_colors_dict(colormap, "rgb")
    elif isinstance(colormap, str) and "rgb" not in colormap and "#" not in colormap:
        if colormap not in clrs.PLOTLY_SCALES.keys():
            raise exceptions.PlotlyError(
                "If 'colormap' is a string, it must be the name "
                "of a Plotly Colorscale. The available colorscale "
                "names are {}".format(clrs.PLOTLY_SCALES.keys())
            )
        else:
            # TODO change below to allow the correct Plotly colorscale
            colormap = clrs.colorscale_to_colors(clrs.PLOTLY_SCALES[colormap])
            # keep only first and last item - fix later
            colormap = [colormap[0]] + [colormap[-1]]
        colormap = clrs.validate_colors(colormap, "rgb")
    else:
        colormap = clrs.validate_colors(colormap, "rgb")

    if not index:
        for name in df:
            headers.append(name)
        for name in headers:
            dataframe.append(df[name].values.tolist())
        # Check for same data-type in df columns
        utils.validate_dataframe(dataframe)
        figure = scatterplot(
            dataframe, headers, diag, size, height, width, title, **kwargs
        )
        return figure
    else:
        # Validate index selection
        if index not in df:
            raise exceptions.PlotlyError(
                "Make sure you set the index "
                "input variable to one of the "
                "column names of your "
                "dataframe."
            )
        index_vals = df[index].values.tolist()
        for name in df:
            if name != index:
                headers.append(name)
        for name in headers:
            dataframe.append(df[name].values.tolist())

        # check for same data-type in each df column
        utils.validate_dataframe(dataframe)
        utils.validate_index(index_vals)

        # check if all colormap keys are in the index
        # if colormap is a dictionary
        if isinstance(colormap, dict):
            for key in colormap:
                if not all(index in colormap for index in index_vals):
                    raise exceptions.PlotlyError(
                        "If colormap is a "
                        "dictionary, all the "
                        "names in the index "
                        "must be keys."
                    )
            figure = scatterplot_dict(
                dataframe,
                headers,
                diag,
                size,
                height,
                width,
                title,
                index,
                index_vals,
                endpts,
                colormap,
                colormap_type,
                **kwargs,
            )
            return figure

        else:
            figure = scatterplot_theme(
                dataframe,
                headers,
                diag,
                size,
                height,
                width,
                title,
                index,
                index_vals,
                endpts,
                colormap,
                colormap_type,
                **kwargs,
            )
            return figure
