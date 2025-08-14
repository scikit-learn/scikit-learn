from plotly import exceptions, optional_imports
from plotly.figure_factory import utils
from plotly.graph_objs import graph_objs

# Optional imports, may be None for users that only use our core functionality.
np = optional_imports.get_module("numpy")
pd = optional_imports.get_module("pandas")
scipy = optional_imports.get_module("scipy")
scipy_stats = optional_imports.get_module("scipy.stats")


DEFAULT_HISTNORM = "probability density"
ALTERNATIVE_HISTNORM = "probability"


def validate_distplot(hist_data, curve_type):
    """
    Distplot-specific validations

    :raises: (PlotlyError) If hist_data is not a list of lists
    :raises: (PlotlyError) If curve_type is not valid (i.e. not 'kde' or
        'normal').
    """
    hist_data_types = (list,)
    if np:
        hist_data_types += (np.ndarray,)
    if pd:
        hist_data_types += (pd.core.series.Series,)

    if not isinstance(hist_data[0], hist_data_types):
        raise exceptions.PlotlyError(
            "Oops, this function was written "
            "to handle multiple datasets, if "
            "you want to plot just one, make "
            "sure your hist_data variable is "
            "still a list of lists, i.e. x = "
            "[1, 2, 3] -> x = [[1, 2, 3]]"
        )

    curve_opts = ("kde", "normal")
    if curve_type not in curve_opts:
        raise exceptions.PlotlyError("curve_type must be defined as 'kde' or 'normal'")

    if not scipy:
        raise ImportError("FigureFactory.create_distplot requires scipy")


def create_distplot(
    hist_data,
    group_labels,
    bin_size=1.0,
    curve_type="kde",
    colors=None,
    rug_text=None,
    histnorm=DEFAULT_HISTNORM,
    show_hist=True,
    show_curve=True,
    show_rug=True,
):
    """
    Function that creates a distplot similar to seaborn.distplot;
    **this function is deprecated**, use instead :mod:`plotly.express`
    functions, for example

    >>> import plotly.express as px
    >>> tips = px.data.tips()
    >>> fig = px.histogram(tips, x="total_bill", y="tip", color="sex", marginal="rug",
    ...                    hover_data=tips.columns)
    >>> fig.show()


    The distplot can be composed of all or any combination of the following
    3 components: (1) histogram, (2) curve: (a) kernel density estimation
    or (b) normal curve, and (3) rug plot. Additionally, multiple distplots
    (from multiple datasets) can be created in the same plot.

    :param (list[list]) hist_data: Use list of lists to plot multiple data
        sets on the same plot.
    :param (list[str]) group_labels: Names for each data set.
    :param (list[float]|float) bin_size: Size of histogram bins.
        Default = 1.
    :param (str) curve_type: 'kde' or 'normal'. Default = 'kde'
    :param (str) histnorm: 'probability density' or 'probability'
        Default = 'probability density'
    :param (bool) show_hist: Add histogram to distplot? Default = True
    :param (bool) show_curve: Add curve to distplot? Default = True
    :param (bool) show_rug: Add rug to distplot? Default = True
    :param (list[str]) colors: Colors for traces.
    :param (list[list]) rug_text: Hovertext values for rug_plot,
    :return (dict): Representation of a distplot figure.

    Example 1: Simple distplot of 1 data set

    >>> from plotly.figure_factory import create_distplot

    >>> hist_data = [[1.1, 1.1, 2.5, 3.0, 3.5,
    ...               3.5, 4.1, 4.4, 4.5, 4.5,
    ...               5.0, 5.0, 5.2, 5.5, 5.5,
    ...               5.5, 5.5, 5.5, 6.1, 7.0]]
    >>> group_labels = ['distplot example']
    >>> fig = create_distplot(hist_data, group_labels)
    >>> fig.show()


    Example 2: Two data sets and added rug text

    >>> from plotly.figure_factory import create_distplot
    >>> # Add histogram data
    >>> hist1_x = [0.8, 1.2, 0.2, 0.6, 1.6,
    ...            -0.9, -0.07, 1.95, 0.9, -0.2,
    ...            -0.5, 0.3, 0.4, -0.37, 0.6]
    >>> hist2_x = [0.8, 1.5, 1.5, 0.6, 0.59,
    ...            1.0, 0.8, 1.7, 0.5, 0.8,
    ...            -0.3, 1.2, 0.56, 0.3, 2.2]

    >>> # Group data together
    >>> hist_data = [hist1_x, hist2_x]

    >>> group_labels = ['2012', '2013']

    >>> # Add text
    >>> rug_text_1 = ['a1', 'b1', 'c1', 'd1', 'e1',
    ...       'f1', 'g1', 'h1', 'i1', 'j1',
    ...       'k1', 'l1', 'm1', 'n1', 'o1']

    >>> rug_text_2 = ['a2', 'b2', 'c2', 'd2', 'e2',
    ...       'f2', 'g2', 'h2', 'i2', 'j2',
    ...       'k2', 'l2', 'm2', 'n2', 'o2']

    >>> # Group text together
    >>> rug_text_all = [rug_text_1, rug_text_2]

    >>> # Create distplot
    >>> fig = create_distplot(
    ...     hist_data, group_labels, rug_text=rug_text_all, bin_size=.2)

    >>> # Add title
    >>> fig.update_layout(title='Dist Plot') # doctest: +SKIP
    >>> fig.show()


    Example 3: Plot with normal curve and hide rug plot

    >>> from plotly.figure_factory import create_distplot
    >>> import numpy as np

    >>> x1 = np.random.randn(190)
    >>> x2 = np.random.randn(200)+1
    >>> x3 = np.random.randn(200)-1
    >>> x4 = np.random.randn(210)+2

    >>> hist_data = [x1, x2, x3, x4]
    >>> group_labels = ['2012', '2013', '2014', '2015']

    >>> fig = create_distplot(
    ...     hist_data, group_labels, curve_type='normal',
    ...     show_rug=False, bin_size=.4)


    Example 4: Distplot with Pandas

    >>> from plotly.figure_factory import create_distplot
    >>> import numpy as np
    >>> import pandas as pd

    >>> df = pd.DataFrame({'2012': np.random.randn(200),
    ...                    '2013': np.random.randn(200)+1})
    >>> fig = create_distplot([df[c] for c in df.columns], df.columns)
    >>> fig.show()
    """
    if colors is None:
        colors = []
    if rug_text is None:
        rug_text = []

    validate_distplot(hist_data, curve_type)
    utils.validate_equal_length(hist_data, group_labels)

    if isinstance(bin_size, (float, int)):
        bin_size = [bin_size] * len(hist_data)

    data = []
    if show_hist:
        hist = _Distplot(
            hist_data,
            histnorm,
            group_labels,
            bin_size,
            curve_type,
            colors,
            rug_text,
            show_hist,
            show_curve,
        ).make_hist()

        data.append(hist)

    if show_curve:
        if curve_type == "normal":
            curve = _Distplot(
                hist_data,
                histnorm,
                group_labels,
                bin_size,
                curve_type,
                colors,
                rug_text,
                show_hist,
                show_curve,
            ).make_normal()
        else:
            curve = _Distplot(
                hist_data,
                histnorm,
                group_labels,
                bin_size,
                curve_type,
                colors,
                rug_text,
                show_hist,
                show_curve,
            ).make_kde()

        data.append(curve)

    if show_rug:
        rug = _Distplot(
            hist_data,
            histnorm,
            group_labels,
            bin_size,
            curve_type,
            colors,
            rug_text,
            show_hist,
            show_curve,
        ).make_rug()

        data.append(rug)
        layout = graph_objs.Layout(
            barmode="overlay",
            hovermode="closest",
            legend=dict(traceorder="reversed"),
            xaxis1=dict(domain=[0.0, 1.0], anchor="y2", zeroline=False),
            yaxis1=dict(domain=[0.35, 1], anchor="free", position=0.0),
            yaxis2=dict(domain=[0, 0.25], anchor="x1", dtick=1, showticklabels=False),
        )
    else:
        layout = graph_objs.Layout(
            barmode="overlay",
            hovermode="closest",
            legend=dict(traceorder="reversed"),
            xaxis1=dict(domain=[0.0, 1.0], anchor="y2", zeroline=False),
            yaxis1=dict(domain=[0.0, 1], anchor="free", position=0.0),
        )

    data = sum(data, [])
    return graph_objs.Figure(data=data, layout=layout)


class _Distplot(object):
    """
    Refer to TraceFactory.create_distplot() for docstring
    """

    def __init__(
        self,
        hist_data,
        histnorm,
        group_labels,
        bin_size,
        curve_type,
        colors,
        rug_text,
        show_hist,
        show_curve,
    ):
        self.hist_data = hist_data
        self.histnorm = histnorm
        self.group_labels = group_labels
        self.bin_size = bin_size
        self.show_hist = show_hist
        self.show_curve = show_curve
        self.trace_number = len(hist_data)
        if rug_text:
            self.rug_text = rug_text
        else:
            self.rug_text = [None] * self.trace_number

        self.start = []
        self.end = []
        if colors:
            self.colors = colors
        else:
            self.colors = [
                "rgb(31, 119, 180)",
                "rgb(255, 127, 14)",
                "rgb(44, 160, 44)",
                "rgb(214, 39, 40)",
                "rgb(148, 103, 189)",
                "rgb(140, 86, 75)",
                "rgb(227, 119, 194)",
                "rgb(127, 127, 127)",
                "rgb(188, 189, 34)",
                "rgb(23, 190, 207)",
            ]
        self.curve_x = [None] * self.trace_number
        self.curve_y = [None] * self.trace_number

        for trace in self.hist_data:
            self.start.append(min(trace) * 1.0)
            self.end.append(max(trace) * 1.0)

    def make_hist(self):
        """
        Makes the histogram(s) for FigureFactory.create_distplot().

        :rtype (list) hist: list of histogram representations
        """
        hist = [None] * self.trace_number

        for index in range(self.trace_number):
            hist[index] = dict(
                type="histogram",
                x=self.hist_data[index],
                xaxis="x1",
                yaxis="y1",
                histnorm=self.histnorm,
                name=self.group_labels[index],
                legendgroup=self.group_labels[index],
                marker=dict(color=self.colors[index % len(self.colors)]),
                autobinx=False,
                xbins=dict(
                    start=self.start[index],
                    end=self.end[index],
                    size=self.bin_size[index],
                ),
                opacity=0.7,
            )
        return hist

    def make_kde(self):
        """
        Makes the kernel density estimation(s) for create_distplot().

        This is called when curve_type = 'kde' in create_distplot().

        :rtype (list) curve: list of kde representations
        """
        curve = [None] * self.trace_number
        for index in range(self.trace_number):
            self.curve_x[index] = [
                self.start[index] + x * (self.end[index] - self.start[index]) / 500
                for x in range(500)
            ]
            self.curve_y[index] = scipy_stats.gaussian_kde(self.hist_data[index])(
                self.curve_x[index]
            )

            if self.histnorm == ALTERNATIVE_HISTNORM:
                self.curve_y[index] *= self.bin_size[index]

        for index in range(self.trace_number):
            curve[index] = dict(
                type="scatter",
                x=self.curve_x[index],
                y=self.curve_y[index],
                xaxis="x1",
                yaxis="y1",
                mode="lines",
                name=self.group_labels[index],
                legendgroup=self.group_labels[index],
                showlegend=False if self.show_hist else True,
                marker=dict(color=self.colors[index % len(self.colors)]),
            )
        return curve

    def make_normal(self):
        """
        Makes the normal curve(s) for create_distplot().

        This is called when curve_type = 'normal' in create_distplot().

        :rtype (list) curve: list of normal curve representations
        """
        curve = [None] * self.trace_number
        mean = [None] * self.trace_number
        sd = [None] * self.trace_number

        for index in range(self.trace_number):
            mean[index], sd[index] = scipy_stats.norm.fit(self.hist_data[index])
            self.curve_x[index] = [
                self.start[index] + x * (self.end[index] - self.start[index]) / 500
                for x in range(500)
            ]
            self.curve_y[index] = scipy_stats.norm.pdf(
                self.curve_x[index], loc=mean[index], scale=sd[index]
            )

            if self.histnorm == ALTERNATIVE_HISTNORM:
                self.curve_y[index] *= self.bin_size[index]

        for index in range(self.trace_number):
            curve[index] = dict(
                type="scatter",
                x=self.curve_x[index],
                y=self.curve_y[index],
                xaxis="x1",
                yaxis="y1",
                mode="lines",
                name=self.group_labels[index],
                legendgroup=self.group_labels[index],
                showlegend=False if self.show_hist else True,
                marker=dict(color=self.colors[index % len(self.colors)]),
            )
        return curve

    def make_rug(self):
        """
        Makes the rug plot(s) for create_distplot().

        :rtype (list) rug: list of rug plot representations
        """
        rug = [None] * self.trace_number
        for index in range(self.trace_number):
            rug[index] = dict(
                type="scatter",
                x=self.hist_data[index],
                y=([self.group_labels[index]] * len(self.hist_data[index])),
                xaxis="x1",
                yaxis="y2",
                mode="markers",
                name=self.group_labels[index],
                legendgroup=self.group_labels[index],
                showlegend=(False if self.show_hist or self.show_curve else True),
                text=self.rug_text[index],
                marker=dict(
                    color=self.colors[index % len(self.colors)], symbol="line-ns-open"
                ),
            )
        return rug
