from collections.abc import Sequence

from plotly import exceptions


def is_sequence(obj):
    return isinstance(obj, Sequence) and not isinstance(obj, str)


def validate_index(index_vals):
    """
    Validates if a list contains all numbers or all strings

    :raises: (PlotlyError) If there are any two items in the list whose
        types differ
    """
    from numbers import Number

    if isinstance(index_vals[0], Number):
        if not all(isinstance(item, Number) for item in index_vals):
            raise exceptions.PlotlyError(
                "Error in indexing column. "
                "Make sure all entries of each "
                "column are all numbers or "
                "all strings."
            )

    elif isinstance(index_vals[0], str):
        if not all(isinstance(item, str) for item in index_vals):
            raise exceptions.PlotlyError(
                "Error in indexing column. "
                "Make sure all entries of each "
                "column are all numbers or "
                "all strings."
            )


def validate_dataframe(array):
    """
    Validates all strings or numbers in each dataframe column

    :raises: (PlotlyError) If there are any two items in any list whose
        types differ
    """
    from numbers import Number

    for vector in array:
        if isinstance(vector[0], Number):
            if not all(isinstance(item, Number) for item in vector):
                raise exceptions.PlotlyError(
                    "Error in dataframe. "
                    "Make sure all entries of "
                    "each column are either "
                    "numbers or strings."
                )
        elif isinstance(vector[0], str):
            if not all(isinstance(item, str) for item in vector):
                raise exceptions.PlotlyError(
                    "Error in dataframe. "
                    "Make sure all entries of "
                    "each column are either "
                    "numbers or strings."
                )


def validate_equal_length(*args):
    """
    Validates that data lists or ndarrays are the same length.

    :raises: (PlotlyError) If any data lists are not the same length.
    """
    length = len(args[0])
    if any(len(lst) != length for lst in args):
        raise exceptions.PlotlyError(
            "Oops! Your data lists or ndarrays should be the same length."
        )


def validate_positive_scalars(**kwargs):
    """
    Validates that all values given in key/val pairs are positive.

    Accepts kwargs to improve Exception messages.

    :raises: (PlotlyError) If any value is < 0 or raises.
    """
    for key, val in kwargs.items():
        try:
            if val <= 0:
                raise ValueError("{} must be > 0, got {}".format(key, val))
        except TypeError:
            raise exceptions.PlotlyError("{} must be a number, got {}".format(key, val))


def flatten(array):
    """
    Uses list comprehension to flatten array

    :param (array): An iterable to flatten
    :raises (PlotlyError): If iterable is not nested.
    :rtype (list): The flattened list.
    """
    try:
        return [item for sublist in array for item in sublist]
    except TypeError:
        raise exceptions.PlotlyError(
            "Your data array could not be "
            "flattened! Make sure your data is "
            "entered as lists or ndarrays!"
        )


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


def annotation_dict_for_label(
    text,
    lane,
    num_of_lanes,
    subplot_spacing,
    row_col="col",
    flipped=True,
    right_side=True,
    text_color="#0f0f0f",
):
    """
    Returns annotation dict for label of n labels of a 1xn or nx1 subplot.

    :param (str) text: the text for a label.
    :param (int) lane: the label number for text. From 1 to n inclusive.
    :param (int) num_of_lanes: the number 'n' of rows or columns in subplot.
    :param (float) subplot_spacing: the value for the horizontal_spacing and
        vertical_spacing params in your plotly.tools.make_subplots() call.
    :param (str) row_col: choose whether labels are placed along rows or
        columns.
    :param (bool) flipped: flips text by 90 degrees. Text is printed
        horizontally if set to True and row_col='row', or if False and
        row_col='col'.
    :param (bool) right_side: only applicable if row_col is set to 'row'.
    :param (str) text_color: color of the text.
    """
    temp = (1 - (num_of_lanes - 1) * subplot_spacing) / (num_of_lanes)
    if not flipped:
        xanchor = "center"
        yanchor = "middle"
        if row_col == "col":
            x = (lane - 1) * (temp + subplot_spacing) + 0.5 * temp
            y = 1.03
            textangle = 0
        elif row_col == "row":
            y = (lane - 1) * (temp + subplot_spacing) + 0.5 * temp
            x = 1.03
            textangle = 90
    else:
        if row_col == "col":
            xanchor = "center"
            yanchor = "bottom"
            x = (lane - 1) * (temp + subplot_spacing) + 0.5 * temp
            y = 1.0
            textangle = 270
        elif row_col == "row":
            yanchor = "middle"
            y = (lane - 1) * (temp + subplot_spacing) + 0.5 * temp
            if right_side:
                x = 1.0
                xanchor = "left"
            else:
                x = -0.01
                xanchor = "right"
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
        text=text,
        font=dict(size=13, color=text_color),
    )
    return annotation_dict


def list_of_options(iterable, conj="and", period=True):
    """
    Returns an English listing of objects seperated by commas ','

    For example, ['foo', 'bar', 'baz'] becomes 'foo, bar and baz'
    if the conjunction 'and' is selected.
    """
    if len(iterable) < 2:
        raise exceptions.PlotlyError(
            "Your list or tuple must contain at least 2 items."
        )
    template = (len(iterable) - 2) * "{}, " + "{} " + conj + " {}" + period * "."
    return template.format(*iterable)
