from plotly import exceptions
from plotly.graph_objs import graph_objs
from plotly.figure_factory import utils


# Default colours for finance charts
_DEFAULT_INCREASING_COLOR = "#3D9970"  # http://clrs.cc
_DEFAULT_DECREASING_COLOR = "#FF4136"


def validate_ohlc(open, high, low, close, direction, **kwargs):
    """
    ohlc and candlestick specific validations

    Specifically, this checks that the high value is the greatest value and
    the low value is the lowest value in each unit.

    See FigureFactory.create_ohlc() or FigureFactory.create_candlestick()
    for params

    :raises: (PlotlyError) If the high value is not the greatest value in
        each unit.
    :raises: (PlotlyError) If the low value is not the lowest value in each
        unit.
    :raises: (PlotlyError) If direction is not 'increasing' or 'decreasing'
    """
    for lst in [open, low, close]:
        for index in range(len(high)):
            if high[index] < lst[index]:
                raise exceptions.PlotlyError(
                    "Oops! Looks like some of "
                    "your high values are less "
                    "the corresponding open, "
                    "low, or close values. "
                    "Double check that your data "
                    "is entered in O-H-L-C order"
                )

    for lst in [open, high, close]:
        for index in range(len(low)):
            if low[index] > lst[index]:
                raise exceptions.PlotlyError(
                    "Oops! Looks like some of "
                    "your low values are greater "
                    "than the corresponding high"
                    ", open, or close values. "
                    "Double check that your data "
                    "is entered in O-H-L-C order"
                )

    direction_opts = ("increasing", "decreasing", "both")
    if direction not in direction_opts:
        raise exceptions.PlotlyError(
            "direction must be defined as " "'increasing', 'decreasing', or " "'both'"
        )


def make_increasing_ohlc(open, high, low, close, dates, **kwargs):
    """
    Makes increasing ohlc sticks

    _make_increasing_ohlc() and _make_decreasing_ohlc separate the
    increasing trace from the decreasing trace so kwargs (such as
    color) can be passed separately to increasing or decreasing traces
    when direction is set to 'increasing' or 'decreasing' in
    FigureFactory.create_candlestick()

    :param (list) open: opening values
    :param (list) high: high values
    :param (list) low: low values
    :param (list) close: closing values
    :param (list) dates: list of datetime objects. Default: None
    :param kwargs: kwargs to be passed to increasing trace via
        plotly.graph_objs.Scatter.

    :rtype (trace) ohlc_incr_data: Scatter trace of all increasing ohlc
        sticks.
    """
    (flat_increase_x, flat_increase_y, text_increase) = _OHLC(
        open, high, low, close, dates
    ).get_increase()

    if "name" in kwargs:
        showlegend = True
    else:
        kwargs.setdefault("name", "Increasing")
        showlegend = False

    kwargs.setdefault("line", dict(color=_DEFAULT_INCREASING_COLOR, width=1))
    kwargs.setdefault("text", text_increase)

    ohlc_incr = dict(
        type="scatter",
        x=flat_increase_x,
        y=flat_increase_y,
        mode="lines",
        showlegend=showlegend,
        **kwargs,
    )
    return ohlc_incr


def make_decreasing_ohlc(open, high, low, close, dates, **kwargs):
    """
    Makes decreasing ohlc sticks

    :param (list) open: opening values
    :param (list) high: high values
    :param (list) low: low values
    :param (list) close: closing values
    :param (list) dates: list of datetime objects. Default: None
    :param kwargs: kwargs to be passed to increasing trace via
        plotly.graph_objs.Scatter.

    :rtype (trace) ohlc_decr_data: Scatter trace of all decreasing ohlc
        sticks.
    """
    (flat_decrease_x, flat_decrease_y, text_decrease) = _OHLC(
        open, high, low, close, dates
    ).get_decrease()

    kwargs.setdefault("line", dict(color=_DEFAULT_DECREASING_COLOR, width=1))
    kwargs.setdefault("text", text_decrease)
    kwargs.setdefault("showlegend", False)
    kwargs.setdefault("name", "Decreasing")

    ohlc_decr = dict(
        type="scatter", x=flat_decrease_x, y=flat_decrease_y, mode="lines", **kwargs
    )
    return ohlc_decr


def create_ohlc(open, high, low, close, dates=None, direction="both", **kwargs):
    """
    **deprecated**, use instead the plotly.graph_objects trace
    :class:`plotly.graph_objects.Ohlc`

    :param (list) open: opening values
    :param (list) high: high values
    :param (list) low: low values
    :param (list) close: closing
    :param (list) dates: list of datetime objects. Default: None
    :param (string) direction: direction can be 'increasing', 'decreasing',
        or 'both'. When the direction is 'increasing', the returned figure
        consists of all units where the close value is greater than the
        corresponding open value, and when the direction is 'decreasing',
        the returned figure consists of all units where the close value is
        less than or equal to the corresponding open value. When the
        direction is 'both', both increasing and decreasing units are
        returned. Default: 'both'
    :param kwargs: kwargs passed through plotly.graph_objs.Scatter.
        These kwargs describe other attributes about the ohlc Scatter trace
        such as the color or the legend name. For more information on valid
        kwargs call help(plotly.graph_objs.Scatter)

    :rtype (dict): returns a representation of an ohlc chart figure.

    Example 1: Simple OHLC chart from a Pandas DataFrame

    >>> from plotly.figure_factory import create_ohlc
    >>> from datetime import datetime

    >>> import pandas as pd
    >>> df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv')
    >>> fig = create_ohlc(df['AAPL.Open'], df['AAPL.High'], df['AAPL.Low'], df['AAPL.Close'], dates=df.index)
    >>> fig.show()
    """
    if dates is not None:
        utils.validate_equal_length(open, high, low, close, dates)
    else:
        utils.validate_equal_length(open, high, low, close)
    validate_ohlc(open, high, low, close, direction, **kwargs)

    if direction == "increasing":
        ohlc_incr = make_increasing_ohlc(open, high, low, close, dates, **kwargs)
        data = [ohlc_incr]
    elif direction == "decreasing":
        ohlc_decr = make_decreasing_ohlc(open, high, low, close, dates, **kwargs)
        data = [ohlc_decr]
    else:
        ohlc_incr = make_increasing_ohlc(open, high, low, close, dates, **kwargs)
        ohlc_decr = make_decreasing_ohlc(open, high, low, close, dates, **kwargs)
        data = [ohlc_incr, ohlc_decr]

    layout = graph_objs.Layout(xaxis=dict(zeroline=False), hovermode="closest")

    return graph_objs.Figure(data=data, layout=layout)


class _OHLC(object):
    """
    Refer to FigureFactory.create_ohlc_increase() for docstring.
    """

    def __init__(self, open, high, low, close, dates, **kwargs):
        self.open = open
        self.high = high
        self.low = low
        self.close = close
        self.empty = [None] * len(open)
        self.dates = dates

        self.all_x = []
        self.all_y = []
        self.increase_x = []
        self.increase_y = []
        self.decrease_x = []
        self.decrease_y = []
        self.get_all_xy()
        self.separate_increase_decrease()

    def get_all_xy(self):
        """
        Zip data to create OHLC shape

        OHLC shape: low to high vertical bar with
        horizontal branches for open and close values.
        If dates were added, the smallest date difference is calculated and
        multiplied by .2 to get the length of the open and close branches.
        If no date data was provided, the x-axis is a list of integers and the
        length of the open and close branches is .2.
        """
        self.all_y = list(
            zip(
                self.open,
                self.open,
                self.high,
                self.low,
                self.close,
                self.close,
                self.empty,
            )
        )
        if self.dates is not None:
            date_dif = []
            for i in range(len(self.dates) - 1):
                date_dif.append(self.dates[i + 1] - self.dates[i])
            date_dif_min = (min(date_dif)) / 5
            self.all_x = [
                [x - date_dif_min, x, x, x, x, x + date_dif_min, None]
                for x in self.dates
            ]
        else:
            self.all_x = [
                [x - 0.2, x, x, x, x, x + 0.2, None] for x in range(len(self.open))
            ]

    def separate_increase_decrease(self):
        """
        Separate data into two groups: increase and decrease

        (1) Increase, where close > open and
        (2) Decrease, where close <= open
        """
        for index in range(len(self.open)):
            if self.close[index] is None:
                pass
            elif self.close[index] > self.open[index]:
                self.increase_x.append(self.all_x[index])
                self.increase_y.append(self.all_y[index])
            else:
                self.decrease_x.append(self.all_x[index])
                self.decrease_y.append(self.all_y[index])

    def get_increase(self):
        """
        Flatten increase data and get increase text

        :rtype (list, list, list): flat_increase_x: x-values for the increasing
            trace, flat_increase_y: y=values for the increasing trace and
            text_increase: hovertext for the increasing trace
        """
        flat_increase_x = utils.flatten(self.increase_x)
        flat_increase_y = utils.flatten(self.increase_y)
        text_increase = ("Open", "Open", "High", "Low", "Close", "Close", "") * (
            len(self.increase_x)
        )

        return flat_increase_x, flat_increase_y, text_increase

    def get_decrease(self):
        """
        Flatten decrease data and get decrease text

        :rtype (list, list, list): flat_decrease_x: x-values for the decreasing
            trace, flat_decrease_y: y=values for the decreasing trace and
            text_decrease: hovertext for the decreasing trace
        """
        flat_decrease_x = utils.flatten(self.decrease_x)
        flat_decrease_y = utils.flatten(self.decrease_y)
        text_decrease = ("Open", "Open", "High", "Low", "Close", "Close", "") * (
            len(self.decrease_x)
        )

        return flat_decrease_x, flat_decrease_y, text_decrease
