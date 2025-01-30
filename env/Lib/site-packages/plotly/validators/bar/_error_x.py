import _plotly_utils.basevalidators


class Error_XValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="error_x", parent_name="bar", **kwargs):
        super(Error_XValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "ErrorX"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            array
                Sets the data corresponding the length of each
                error bar. Values are plotted relative to the
                underlying data.
            arrayminus
                Sets the data corresponding the length of each
                error bar in the bottom (left) direction for
                vertical (horizontal) bars Values are plotted
                relative to the underlying data.
            arrayminussrc
                Sets the source reference on Chart Studio Cloud
                for `arrayminus`.
            arraysrc
                Sets the source reference on Chart Studio Cloud
                for `array`.
            color
                Sets the stoke color of the error bars.
            copy_ystyle

            symmetric
                Determines whether or not the error bars have
                the same length in both direction (top/bottom
                for vertical bars, left/right for horizontal
                bars.
            thickness
                Sets the thickness (in px) of the error bars.
            traceref

            tracerefminus

            type
                Determines the rule used to generate the error
                bars. If *constant`, the bar lengths are of a
                constant value. Set this constant in `value`.
                If "percent", the bar lengths correspond to a
                percentage of underlying data. Set this
                percentage in `value`. If "sqrt", the bar
                lengths correspond to the square of the
                underlying data. If "data", the bar lengths are
                set with data set `array`.
            value
                Sets the value of either the percentage (if
                `type` is set to "percent") or the constant (if
                `type` is set to "constant") corresponding to
                the lengths of the error bars.
            valueminus
                Sets the value of either the percentage (if
                `type` is set to "percent") or the constant (if
                `type` is set to "constant") corresponding to
                the lengths of the error bars in the bottom
                (left) direction for vertical (horizontal) bars
            visible
                Determines whether or not this set of error
                bars is visible.
            width
                Sets the width (in px) of the cross-bar at both
                ends of the error bars.
""",
            ),
            **kwargs,
        )
