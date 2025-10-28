from plotly.figure_factory import utils
from plotly.figure_factory._ohlc import (
    _DEFAULT_INCREASING_COLOR,
    _DEFAULT_DECREASING_COLOR,
    validate_ohlc,
)
from plotly.graph_objs import graph_objs


def make_increasing_candle(open, high, low, close, dates, **kwargs):
    """
    Makes boxplot trace for increasing candlesticks

    _make_increasing_candle() and _make_decreasing_candle separate the
    increasing traces from the decreasing traces so kwargs (such as
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

    :rtype (list) candle_incr_data: list of the box trace for
        increasing candlesticks.
    """
    increase_x, increase_y = _Candlestick(
        open, high, low, close, dates, **kwargs
    ).get_candle_increase()

    if "line" in kwargs:
        kwargs.setdefault("fillcolor", kwargs["line"]["color"])
    else:
        kwargs.setdefault("fillcolor", _DEFAULT_INCREASING_COLOR)
    if "name" in kwargs:
        kwargs.setdefault("showlegend", True)
    else:
        kwargs.setdefault("showlegend", False)
    kwargs.setdefault("name", "Increasing")
    kwargs.setdefault("line", dict(color=_DEFAULT_INCREASING_COLOR))

    candle_incr_data = dict(
        type="box",
        x=increase_x,
        y=increase_y,
        whiskerwidth=0,
        boxpoints=False,
        **kwargs,
    )

    return [candle_incr_data]


def make_decreasing_candle(open, high, low, close, dates, **kwargs):
    """
    Makes boxplot trace for decreasing candlesticks

    :param (list) open: opening values
    :param (list) high: high values
    :param (list) low: low values
    :param (list) close: closing values
    :param (list) dates: list of datetime objects. Default: None
    :param kwargs: kwargs to be passed to decreasing trace via
        plotly.graph_objs.Scatter.

    :rtype (list) candle_decr_data: list of the box trace for
        decreasing candlesticks.
    """

    decrease_x, decrease_y = _Candlestick(
        open, high, low, close, dates, **kwargs
    ).get_candle_decrease()

    if "line" in kwargs:
        kwargs.setdefault("fillcolor", kwargs["line"]["color"])
    else:
        kwargs.setdefault("fillcolor", _DEFAULT_DECREASING_COLOR)
    kwargs.setdefault("showlegend", False)
    kwargs.setdefault("line", dict(color=_DEFAULT_DECREASING_COLOR))
    kwargs.setdefault("name", "Decreasing")

    candle_decr_data = dict(
        type="box",
        x=decrease_x,
        y=decrease_y,
        whiskerwidth=0,
        boxpoints=False,
        **kwargs,
    )

    return [candle_decr_data]


def create_candlestick(open, high, low, close, dates=None, direction="both", **kwargs):
    """
    **deprecated**, use instead the plotly.graph_objects trace
    :class:`plotly.graph_objects.Candlestick`

    :param (list) open: opening values
    :param (list) high: high values
    :param (list) low: low values
    :param (list) close: closing values
    :param (list) dates: list of datetime objects. Default: None
    :param (string) direction: direction can be 'increasing', 'decreasing',
        or 'both'. When the direction is 'increasing', the returned figure
        consists of all candlesticks where the close value is greater than
        the corresponding open value, and when the direction is
        'decreasing', the returned figure consists of all candlesticks
        where the close value is less than or equal to the corresponding
        open value. When the direction is 'both', both increasing and
        decreasing candlesticks are returned. Default: 'both'
    :param kwargs: kwargs passed through plotly.graph_objs.Scatter.
        These kwargs describe other attributes about the ohlc Scatter trace
        such as the color or the legend name. For more information on valid
        kwargs call help(plotly.graph_objs.Scatter)

    :rtype (dict): returns a representation of candlestick chart figure.

    Example 1: Simple candlestick chart from a Pandas DataFrame

    >>> from plotly.figure_factory import create_candlestick
    >>> from datetime import datetime
    >>> import pandas as pd

    >>> df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv')
    >>> fig = create_candlestick(df['AAPL.Open'], df['AAPL.High'], df['AAPL.Low'], df['AAPL.Close'],
    ...                          dates=df.index)
    >>> fig.show()

    Example 2: Customize the candlestick colors

    >>> from plotly.figure_factory import create_candlestick
    >>> from plotly.graph_objs import Line, Marker
    >>> from datetime import datetime

    >>> import pandas as pd
    >>> df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv')

    >>> # Make increasing candlesticks and customize their color and name
    >>> fig_increasing = create_candlestick(df['AAPL.Open'], df['AAPL.High'], df['AAPL.Low'], df['AAPL.Close'],
    ...     dates=df.index,
    ...     direction='increasing', name='AAPL',
    ...     marker=Marker(color='rgb(150, 200, 250)'),
    ...     line=Line(color='rgb(150, 200, 250)'))

    >>> # Make decreasing candlesticks and customize their color and name
    >>> fig_decreasing = create_candlestick(df['AAPL.Open'], df['AAPL.High'], df['AAPL.Low'], df['AAPL.Close'],
    ...     dates=df.index,
    ...     direction='decreasing',
    ...     marker=Marker(color='rgb(128, 128, 128)'),
    ...     line=Line(color='rgb(128, 128, 128)'))

    >>> # Initialize the figure
    >>> fig = fig_increasing

    >>> # Add decreasing data with .extend()
    >>> fig.add_trace(fig_decreasing['data']) # doctest: +SKIP
    >>> fig.show()

    Example 3: Candlestick chart with datetime objects

    >>> from plotly.figure_factory import create_candlestick

    >>> from datetime import datetime

    >>> # Add data
    >>> open_data = [33.0, 33.3, 33.5, 33.0, 34.1]
    >>> high_data = [33.1, 33.3, 33.6, 33.2, 34.8]
    >>> low_data = [32.7, 32.7, 32.8, 32.6, 32.8]
    >>> close_data = [33.0, 32.9, 33.3, 33.1, 33.1]
    >>> dates = [datetime(year=2013, month=10, day=10),
    ...          datetime(year=2013, month=11, day=10),
    ...          datetime(year=2013, month=12, day=10),
    ...          datetime(year=2014, month=1, day=10),
    ...          datetime(year=2014, month=2, day=10)]

    >>> # Create ohlc
    >>> fig = create_candlestick(open_data, high_data,
    ...     low_data, close_data, dates=dates)
    >>> fig.show()
    """
    if dates is not None:
        utils.validate_equal_length(open, high, low, close, dates)
    else:
        utils.validate_equal_length(open, high, low, close)
    validate_ohlc(open, high, low, close, direction, **kwargs)

    if direction == "increasing":
        candle_incr_data = make_increasing_candle(
            open, high, low, close, dates, **kwargs
        )
        data = candle_incr_data
    elif direction == "decreasing":
        candle_decr_data = make_decreasing_candle(
            open, high, low, close, dates, **kwargs
        )
        data = candle_decr_data
    else:
        candle_incr_data = make_increasing_candle(
            open, high, low, close, dates, **kwargs
        )
        candle_decr_data = make_decreasing_candle(
            open, high, low, close, dates, **kwargs
        )
        data = candle_incr_data + candle_decr_data

    layout = graph_objs.Layout()
    return graph_objs.Figure(data=data, layout=layout)


class _Candlestick(object):
    """
    Refer to FigureFactory.create_candlestick() for docstring.
    """

    def __init__(self, open, high, low, close, dates, **kwargs):
        self.open = open
        self.high = high
        self.low = low
        self.close = close
        if dates is not None:
            self.x = dates
        else:
            self.x = [x for x in range(len(self.open))]
        self.get_candle_increase()

    def get_candle_increase(self):
        """
        Separate increasing data from decreasing data.

        The data is increasing when close value > open value
        and decreasing when the close value <= open value.
        """
        increase_y = []
        increase_x = []
        for index in range(len(self.open)):
            if self.close[index] > self.open[index]:
                increase_y.append(self.low[index])
                increase_y.append(self.open[index])
                increase_y.append(self.close[index])
                increase_y.append(self.close[index])
                increase_y.append(self.close[index])
                increase_y.append(self.high[index])
                increase_x.append(self.x[index])

        increase_x = [[x, x, x, x, x, x] for x in increase_x]
        increase_x = utils.flatten(increase_x)

        return increase_x, increase_y

    def get_candle_decrease(self):
        """
        Separate increasing data from decreasing data.

        The data is increasing when close value > open value
        and decreasing when the close value <= open value.
        """
        decrease_y = []
        decrease_x = []
        for index in range(len(self.open)):
            if self.close[index] <= self.open[index]:
                decrease_y.append(self.low[index])
                decrease_y.append(self.open[index])
                decrease_y.append(self.close[index])
                decrease_y.append(self.close[index])
                decrease_y.append(self.close[index])
                decrease_y.append(self.high[index])
                decrease_x.append(self.x[index])

        decrease_x = [[x, x, x, x, x, x] for x in decrease_x]
        decrease_x = utils.flatten(decrease_x)

        return decrease_x, decrease_y
