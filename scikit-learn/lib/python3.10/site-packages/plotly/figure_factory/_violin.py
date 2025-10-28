from numbers import Number

from plotly import exceptions, optional_imports
import plotly.colors as clrs
from plotly.graph_objs import graph_objs
from plotly.subplots import make_subplots

pd = optional_imports.get_module("pandas")
np = optional_imports.get_module("numpy")
scipy_stats = optional_imports.get_module("scipy.stats")


def calc_stats(data):
    """
    Calculate statistics for use in violin plot.
    """
    x = np.asarray(data, float)
    vals_min = np.min(x)
    vals_max = np.max(x)
    q2 = np.percentile(x, 50, interpolation="linear")
    q1 = np.percentile(x, 25, interpolation="lower")
    q3 = np.percentile(x, 75, interpolation="higher")
    iqr = q3 - q1
    whisker_dist = 1.5 * iqr

    # in order to prevent drawing whiskers outside the interval
    # of data one defines the whisker positions as:
    d1 = np.min(x[x >= (q1 - whisker_dist)])
    d2 = np.max(x[x <= (q3 + whisker_dist)])
    return {
        "min": vals_min,
        "max": vals_max,
        "q1": q1,
        "q2": q2,
        "q3": q3,
        "d1": d1,
        "d2": d2,
    }


def make_half_violin(x, y, fillcolor="#1f77b4", linecolor="rgb(0, 0, 0)"):
    """
    Produces a sideways probability distribution fig violin plot.
    """
    text = [
        "(pdf(y), y)=(" + "{:0.2f}".format(x[i]) + ", " + "{:0.2f}".format(y[i]) + ")"
        for i in range(len(x))
    ]

    return graph_objs.Scatter(
        x=x,
        y=y,
        mode="lines",
        name="",
        text=text,
        fill="tonextx",
        fillcolor=fillcolor,
        line=graph_objs.scatter.Line(width=0.5, color=linecolor, shape="spline"),
        hoverinfo="text",
        opacity=0.5,
    )


def make_violin_rugplot(vals, pdf_max, distance, color="#1f77b4"):
    """
    Returns a rugplot fig for a violin plot.
    """
    return graph_objs.Scatter(
        y=vals,
        x=[-pdf_max - distance] * len(vals),
        marker=graph_objs.scatter.Marker(color=color, symbol="line-ew-open"),
        mode="markers",
        name="",
        showlegend=False,
        hoverinfo="y",
    )


def make_non_outlier_interval(d1, d2):
    """
    Returns the scatterplot fig of most of a violin plot.
    """
    return graph_objs.Scatter(
        x=[0, 0],
        y=[d1, d2],
        name="",
        mode="lines",
        line=graph_objs.scatter.Line(width=1.5, color="rgb(0,0,0)"),
    )


def make_quartiles(q1, q3):
    """
    Makes the upper and lower quartiles for a violin plot.
    """
    return graph_objs.Scatter(
        x=[0, 0],
        y=[q1, q3],
        text=[
            "lower-quartile: " + "{:0.2f}".format(q1),
            "upper-quartile: " + "{:0.2f}".format(q3),
        ],
        mode="lines",
        line=graph_objs.scatter.Line(width=4, color="rgb(0,0,0)"),
        hoverinfo="text",
    )


def make_median(q2):
    """
    Formats the 'median' hovertext for a violin plot.
    """
    return graph_objs.Scatter(
        x=[0],
        y=[q2],
        text=["median: " + "{:0.2f}".format(q2)],
        mode="markers",
        marker=dict(symbol="square", color="rgb(255,255,255)"),
        hoverinfo="text",
    )


def make_XAxis(xaxis_title, xaxis_range):
    """
    Makes the x-axis for a violin plot.
    """
    xaxis = graph_objs.layout.XAxis(
        title=xaxis_title,
        range=xaxis_range,
        showgrid=False,
        zeroline=False,
        showline=False,
        mirror=False,
        ticks="",
        showticklabels=False,
    )
    return xaxis


def make_YAxis(yaxis_title):
    """
    Makes the y-axis for a violin plot.
    """
    yaxis = graph_objs.layout.YAxis(
        title=yaxis_title,
        showticklabels=True,
        autorange=True,
        ticklen=4,
        showline=True,
        zeroline=False,
        showgrid=False,
        mirror=False,
    )
    return yaxis


def violinplot(vals, fillcolor="#1f77b4", rugplot=True):
    """
    Refer to FigureFactory.create_violin() for docstring.
    """
    vals = np.asarray(vals, float)
    #  summary statistics
    vals_min = calc_stats(vals)["min"]
    vals_max = calc_stats(vals)["max"]
    q1 = calc_stats(vals)["q1"]
    q2 = calc_stats(vals)["q2"]
    q3 = calc_stats(vals)["q3"]
    d1 = calc_stats(vals)["d1"]
    d2 = calc_stats(vals)["d2"]

    # kernel density estimation of pdf
    pdf = scipy_stats.gaussian_kde(vals)
    # grid over the data interval
    xx = np.linspace(vals_min, vals_max, 100)
    # evaluate the pdf at the grid xx
    yy = pdf(xx)
    max_pdf = np.max(yy)
    # distance from the violin plot to rugplot
    distance = (2.0 * max_pdf) / 10 if rugplot else 0
    # range for x values in the plot
    plot_xrange = [-max_pdf - distance - 0.1, max_pdf + 0.1]
    plot_data = [
        make_half_violin(-yy, xx, fillcolor=fillcolor),
        make_half_violin(yy, xx, fillcolor=fillcolor),
        make_non_outlier_interval(d1, d2),
        make_quartiles(q1, q3),
        make_median(q2),
    ]
    if rugplot:
        plot_data.append(
            make_violin_rugplot(vals, max_pdf, distance=distance, color=fillcolor)
        )
    return plot_data, plot_xrange


def violin_no_colorscale(
    data,
    data_header,
    group_header,
    colors,
    use_colorscale,
    group_stats,
    rugplot,
    sort,
    height,
    width,
    title,
):
    """
    Refer to FigureFactory.create_violin() for docstring.

    Returns fig for violin plot without colorscale.

    """

    # collect all group names
    group_name = []
    for name in data[group_header]:
        if name not in group_name:
            group_name.append(name)
    if sort:
        group_name.sort()

    gb = data.groupby([group_header])
    L = len(group_name)

    fig = make_subplots(
        rows=1, cols=L, shared_yaxes=True, horizontal_spacing=0.025, print_grid=False
    )
    color_index = 0
    for k, gr in enumerate(group_name):
        vals = np.asarray(gb.get_group(gr)[data_header], float)
        if color_index >= len(colors):
            color_index = 0
        plot_data, plot_xrange = violinplot(
            vals, fillcolor=colors[color_index], rugplot=rugplot
        )
        for item in plot_data:
            fig.append_trace(item, 1, k + 1)
        color_index += 1

        # add violin plot labels
        fig["layout"].update(
            {"xaxis{}".format(k + 1): make_XAxis(group_name[k], plot_xrange)}
        )

    # set the sharey axis style
    fig["layout"].update({"yaxis{}".format(1): make_YAxis("")})
    fig["layout"].update(
        title=title,
        showlegend=False,
        hovermode="closest",
        autosize=False,
        height=height,
        width=width,
    )

    return fig


def violin_colorscale(
    data,
    data_header,
    group_header,
    colors,
    use_colorscale,
    group_stats,
    rugplot,
    sort,
    height,
    width,
    title,
):
    """
    Refer to FigureFactory.create_violin() for docstring.

    Returns fig for violin plot with colorscale.

    """

    # collect all group names
    group_name = []
    for name in data[group_header]:
        if name not in group_name:
            group_name.append(name)
    if sort:
        group_name.sort()

    # make sure all group names are keys in group_stats
    for group in group_name:
        if group not in group_stats:
            raise exceptions.PlotlyError(
                "All values/groups in the index "
                "column must be represented "
                "as a key in group_stats."
            )

    gb = data.groupby([group_header])
    L = len(group_name)

    fig = make_subplots(
        rows=1, cols=L, shared_yaxes=True, horizontal_spacing=0.025, print_grid=False
    )

    # prepare low and high color for colorscale
    lowcolor = clrs.color_parser(colors[0], clrs.unlabel_rgb)
    highcolor = clrs.color_parser(colors[1], clrs.unlabel_rgb)

    # find min and max values in group_stats
    group_stats_values = []
    for key in group_stats:
        group_stats_values.append(group_stats[key])

    max_value = max(group_stats_values)
    min_value = min(group_stats_values)

    for k, gr in enumerate(group_name):
        vals = np.asarray(gb.get_group(gr)[data_header], float)

        # find intermediate color from colorscale
        intermed = (group_stats[gr] - min_value) / (max_value - min_value)
        intermed_color = clrs.find_intermediate_color(lowcolor, highcolor, intermed)

        plot_data, plot_xrange = violinplot(
            vals, fillcolor="rgb{}".format(intermed_color), rugplot=rugplot
        )
        for item in plot_data:
            fig.append_trace(item, 1, k + 1)
        fig["layout"].update(
            {"xaxis{}".format(k + 1): make_XAxis(group_name[k], plot_xrange)}
        )
    # add colorbar to plot
    trace_dummy = graph_objs.Scatter(
        x=[0],
        y=[0],
        mode="markers",
        marker=dict(
            size=2,
            cmin=min_value,
            cmax=max_value,
            colorscale=[[0, colors[0]], [1, colors[1]]],
            showscale=True,
        ),
        showlegend=False,
    )
    fig.append_trace(trace_dummy, 1, L)

    # set the sharey axis style
    fig["layout"].update({"yaxis{}".format(1): make_YAxis("")})
    fig["layout"].update(
        title=title,
        showlegend=False,
        hovermode="closest",
        autosize=False,
        height=height,
        width=width,
    )

    return fig


def violin_dict(
    data,
    data_header,
    group_header,
    colors,
    use_colorscale,
    group_stats,
    rugplot,
    sort,
    height,
    width,
    title,
):
    """
    Refer to FigureFactory.create_violin() for docstring.

    Returns fig for violin plot without colorscale.

    """

    # collect all group names
    group_name = []
    for name in data[group_header]:
        if name not in group_name:
            group_name.append(name)

    if sort:
        group_name.sort()

    # check if all group names appear in colors dict
    for group in group_name:
        if group not in colors:
            raise exceptions.PlotlyError(
                "If colors is a dictionary, all "
                "the group names must appear as "
                "keys in colors."
            )

    gb = data.groupby([group_header])
    L = len(group_name)

    fig = make_subplots(
        rows=1, cols=L, shared_yaxes=True, horizontal_spacing=0.025, print_grid=False
    )

    for k, gr in enumerate(group_name):
        vals = np.asarray(gb.get_group(gr)[data_header], float)
        plot_data, plot_xrange = violinplot(vals, fillcolor=colors[gr], rugplot=rugplot)
        for item in plot_data:
            fig.append_trace(item, 1, k + 1)

        # add violin plot labels
        fig["layout"].update(
            {"xaxis{}".format(k + 1): make_XAxis(group_name[k], plot_xrange)}
        )

    # set the sharey axis style
    fig["layout"].update({"yaxis{}".format(1): make_YAxis("")})
    fig["layout"].update(
        title=title,
        showlegend=False,
        hovermode="closest",
        autosize=False,
        height=height,
        width=width,
    )

    return fig


def create_violin(
    data,
    data_header=None,
    group_header=None,
    colors=None,
    use_colorscale=False,
    group_stats=None,
    rugplot=True,
    sort=False,
    height=450,
    width=600,
    title="Violin and Rug Plot",
):
    """
    **deprecated**, use instead the plotly.graph_objects trace
    :class:`plotly.graph_objects.Violin`.

    :param (list|array) data: accepts either a list of numerical values,
        a list of dictionaries all with identical keys and at least one
        column of numeric values, or a pandas dataframe with at least one
        column of numbers.
    :param (str) data_header: the header of the data column to be used
        from an inputted pandas dataframe. Not applicable if 'data' is
        a list of numeric values.
    :param (str) group_header: applicable if grouping data by a variable.
        'group_header' must be set to the name of the grouping variable.
    :param (str|tuple|list|dict) colors: either a plotly scale name,
        an rgb or hex color, a color tuple, a list of colors or a
        dictionary. An rgb color is of the form 'rgb(x, y, z)' where
        x, y and z belong to the interval [0, 255] and a color tuple is a
        tuple of the form (a, b, c) where a, b and c belong to [0, 1].
        If colors is a list, it must contain valid color types as its
        members.
    :param (bool) use_colorscale: only applicable if grouping by another
        variable. Will implement a colorscale based on the first 2 colors
        of param colors. This means colors must be a list with at least 2
        colors in it (Plotly colorscales are accepted since they map to a
        list of two rgb colors). Default = False
    :param (dict) group_stats: a dictionary where each key is a unique
        value from the group_header column in data. Each value must be a
        number and will be used to color the violin plots if a colorscale
        is being used.
    :param (bool) rugplot: determines if a rugplot is draw on violin plot.
        Default = True
    :param (bool) sort: determines if violins are sorted
        alphabetically (True) or by input order (False). Default = False
    :param (float) height: the height of the violin plot.
    :param (float) width: the width of the violin plot.
    :param (str) title: the title of the violin plot.

    Example 1: Single Violin Plot

    >>> from plotly.figure_factory import create_violin
    >>> import plotly.graph_objs as graph_objects

    >>> import numpy as np
    >>> from scipy import stats

    >>> # create list of random values
    >>> data_list = np.random.randn(100)

    >>> # create violin fig
    >>> fig = create_violin(data_list, colors='#604d9e')

    >>> # plot
    >>> fig.show()

    Example 2: Multiple Violin Plots with Qualitative Coloring

    >>> from plotly.figure_factory import create_violin
    >>> import plotly.graph_objs as graph_objects

    >>> import numpy as np
    >>> import pandas as pd
    >>> from scipy import stats

    >>> # create dataframe
    >>> np.random.seed(619517)
    >>> Nr=250
    >>> y = np.random.randn(Nr)
    >>> gr = np.random.choice(list("ABCDE"), Nr)
    >>> norm_params=[(0, 1.2), (0.7, 1), (-0.5, 1.4), (0.3, 1), (0.8, 0.9)]

    >>> for i, letter in enumerate("ABCDE"):
    ...     y[gr == letter] *=norm_params[i][1]+ norm_params[i][0]
    >>> df = pd.DataFrame(dict(Score=y, Group=gr))

    >>> # create violin fig
    >>> fig = create_violin(df, data_header='Score', group_header='Group',
    ...                    sort=True, height=600, width=1000)

    >>> # plot
    >>> fig.show()

    Example 3: Violin Plots with Colorscale

    >>> from plotly.figure_factory import create_violin
    >>> import plotly.graph_objs as graph_objects

    >>> import numpy as np
    >>> import pandas as pd
    >>> from scipy import stats

    >>> # create dataframe
    >>> np.random.seed(619517)
    >>> Nr=250
    >>> y = np.random.randn(Nr)
    >>> gr = np.random.choice(list("ABCDE"), Nr)
    >>> norm_params=[(0, 1.2), (0.7, 1), (-0.5, 1.4), (0.3, 1), (0.8, 0.9)]

    >>> for i, letter in enumerate("ABCDE"):
    ...     y[gr == letter] *=norm_params[i][1]+ norm_params[i][0]
    >>> df = pd.DataFrame(dict(Score=y, Group=gr))

    >>> # define header params
    >>> data_header = 'Score'
    >>> group_header = 'Group'

    >>> # make groupby object with pandas
    >>> group_stats = {}
    >>> groupby_data = df.groupby([group_header])

    >>> for group in "ABCDE":
    ...     data_from_group = groupby_data.get_group(group)[data_header]
    ...     # take a stat of the grouped data
    ...     stat = np.median(data_from_group)
    ...     # add to dictionary
    ...     group_stats[group] = stat

    >>> # create violin fig
    >>> fig = create_violin(df, data_header='Score', group_header='Group',
    ...                     height=600, width=1000, use_colorscale=True,
    ...                     group_stats=group_stats)

    >>> # plot
    >>> fig.show()
    """

    # Validate colors
    if isinstance(colors, dict):
        valid_colors = clrs.validate_colors_dict(colors, "rgb")
    else:
        valid_colors = clrs.validate_colors(colors, "rgb")

    # validate data and choose plot type
    if group_header is None:
        if isinstance(data, list):
            if len(data) <= 0:
                raise exceptions.PlotlyError(
                    "If data is a list, it must be "
                    "nonempty and contain either "
                    "numbers or dictionaries."
                )

            if not all(isinstance(element, Number) for element in data):
                raise exceptions.PlotlyError(
                    "If data is a list, it must contain only numbers."
                )

        if pd and isinstance(data, pd.core.frame.DataFrame):
            if data_header is None:
                raise exceptions.PlotlyError(
                    "data_header must be the "
                    "column name with the "
                    "desired numeric data for "
                    "the violin plot."
                )

            data = data[data_header].values.tolist()

        # call the plotting functions
        plot_data, plot_xrange = violinplot(
            data, fillcolor=valid_colors[0], rugplot=rugplot
        )

        layout = graph_objs.Layout(
            title=title,
            autosize=False,
            font=graph_objs.layout.Font(size=11),
            height=height,
            showlegend=False,
            width=width,
            xaxis=make_XAxis("", plot_xrange),
            yaxis=make_YAxis(""),
            hovermode="closest",
        )
        layout["yaxis"].update(dict(showline=False, showticklabels=False, ticks=""))

        fig = graph_objs.Figure(data=plot_data, layout=layout)

        return fig

    else:
        if not isinstance(data, pd.core.frame.DataFrame):
            raise exceptions.PlotlyError(
                "Error. You must use a pandas "
                "DataFrame if you are using a "
                "group header."
            )

        if data_header is None:
            raise exceptions.PlotlyError(
                "data_header must be the column "
                "name with the desired numeric "
                "data for the violin plot."
            )

        if use_colorscale is False:
            if isinstance(valid_colors, dict):
                # validate colors dict choice below
                fig = violin_dict(
                    data,
                    data_header,
                    group_header,
                    valid_colors,
                    use_colorscale,
                    group_stats,
                    rugplot,
                    sort,
                    height,
                    width,
                    title,
                )
                return fig
            else:
                fig = violin_no_colorscale(
                    data,
                    data_header,
                    group_header,
                    valid_colors,
                    use_colorscale,
                    group_stats,
                    rugplot,
                    sort,
                    height,
                    width,
                    title,
                )
                return fig
        else:
            if isinstance(valid_colors, dict):
                raise exceptions.PlotlyError(
                    "The colors param cannot be "
                    "a dictionary if you are "
                    "using a colorscale."
                )

            if len(valid_colors) < 2:
                raise exceptions.PlotlyError(
                    "colors must be a list with "
                    "at least 2 colors. A "
                    "Plotly scale is allowed."
                )

            if not isinstance(group_stats, dict):
                raise exceptions.PlotlyError(
                    "Your group_stats param must be a dictionary."
                )

            fig = violin_colorscale(
                data,
                data_header,
                group_header,
                valid_colors,
                use_colorscale,
                group_stats,
                rugplot,
                sort,
                height,
                width,
                title,
            )
            return fig
