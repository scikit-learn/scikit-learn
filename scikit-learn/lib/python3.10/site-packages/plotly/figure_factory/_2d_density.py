from numbers import Number

import plotly.exceptions

import plotly.colors as clrs
from plotly.graph_objs import graph_objs


def make_linear_colorscale(colors):
    """
    Makes a list of colors into a colorscale-acceptable form

    For documentation regarding to the form of the output, see
    https://plot.ly/python/reference/#mesh3d-colorscale
    """
    scale = 1.0 / (len(colors) - 1)
    return [[i * scale, color] for i, color in enumerate(colors)]


def create_2d_density(
    x,
    y,
    colorscale="Earth",
    ncontours=20,
    hist_color=(0, 0, 0.5),
    point_color=(0, 0, 0.5),
    point_size=2,
    title="2D Density Plot",
    height=600,
    width=600,
):
    """
    **deprecated**, use instead
    :func:`plotly.express.density_heatmap`.

    :param (list|array) x: x-axis data for plot generation
    :param (list|array) y: y-axis data for plot generation
    :param (str|tuple|list) colorscale: either a plotly scale name, an rgb
        or hex color, a color tuple or a list or tuple of colors. An rgb
        color is of the form 'rgb(x, y, z)' where x, y, z belong to the
        interval [0, 255] and a color tuple is a tuple of the form
        (a, b, c) where a, b and c belong to [0, 1]. If colormap is a
        list, it must contain the valid color types aforementioned as its
        members.
    :param (int) ncontours: the number of 2D contours to draw on the plot
    :param (str) hist_color: the color of the plotted histograms
    :param (str) point_color: the color of the scatter points
    :param (str) point_size: the color of the scatter points
    :param (str) title: set the title for the plot
    :param (float) height: the height of the chart
    :param (float) width: the width of the chart

    Examples
    --------

    Example 1: Simple 2D Density Plot

    >>> from plotly.figure_factory import create_2d_density
    >>> import numpy as np

    >>> # Make data points
    >>> t = np.linspace(-1,1.2,2000)
    >>> x = (t**3)+(0.3*np.random.randn(2000))
    >>> y = (t**6)+(0.3*np.random.randn(2000))

    >>> # Create a figure
    >>> fig = create_2d_density(x, y)

    >>> # Plot the data
    >>> fig.show()

    Example 2: Using Parameters

    >>> from plotly.figure_factory import create_2d_density

    >>> import numpy as np

    >>> # Make data points
    >>> t = np.linspace(-1,1.2,2000)
    >>> x = (t**3)+(0.3*np.random.randn(2000))
    >>> y = (t**6)+(0.3*np.random.randn(2000))

    >>> # Create custom colorscale
    >>> colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)',
    ...              (1, 1, 0.2), (0.98,0.98,0.98)]

    >>> # Create a figure
    >>> fig = create_2d_density(x, y, colorscale=colorscale,
    ...       hist_color='rgb(255, 237, 222)', point_size=3)

    >>> # Plot the data
    >>> fig.show()
    """

    # validate x and y are filled with numbers only
    for array in [x, y]:
        if not all(isinstance(element, Number) for element in array):
            raise plotly.exceptions.PlotlyError(
                "All elements of your 'x' and 'y' lists must be numbers."
            )

    # validate x and y are the same length
    if len(x) != len(y):
        raise plotly.exceptions.PlotlyError(
            "Both lists 'x' and 'y' must be the same length."
        )

    colorscale = clrs.validate_colors(colorscale, "rgb")
    colorscale = make_linear_colorscale(colorscale)

    # validate hist_color and point_color
    hist_color = clrs.validate_colors(hist_color, "rgb")
    point_color = clrs.validate_colors(point_color, "rgb")

    trace1 = graph_objs.Scatter(
        x=x,
        y=y,
        mode="markers",
        name="points",
        marker=dict(color=point_color[0], size=point_size, opacity=0.4),
    )
    trace2 = graph_objs.Histogram2dContour(
        x=x,
        y=y,
        name="density",
        ncontours=ncontours,
        colorscale=colorscale,
        reversescale=True,
        showscale=False,
    )
    trace3 = graph_objs.Histogram(
        x=x, name="x density", marker=dict(color=hist_color[0]), yaxis="y2"
    )
    trace4 = graph_objs.Histogram(
        y=y, name="y density", marker=dict(color=hist_color[0]), xaxis="x2"
    )
    data = [trace1, trace2, trace3, trace4]

    layout = graph_objs.Layout(
        showlegend=False,
        autosize=False,
        title=title,
        height=height,
        width=width,
        xaxis=dict(domain=[0, 0.85], showgrid=False, zeroline=False),
        yaxis=dict(domain=[0, 0.85], showgrid=False, zeroline=False),
        margin=dict(t=50),
        hovermode="closest",
        bargap=0,
        xaxis2=dict(domain=[0.85, 1], showgrid=False, zeroline=False),
        yaxis2=dict(domain=[0.85, 1], showgrid=False, zeroline=False),
    )

    fig = graph_objs.Figure(data=data, layout=layout)
    return fig
