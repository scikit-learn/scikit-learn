import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="line", parent_name="choroplethmapbox.marker", **kwargs
    ):
        super(LineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Line"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the marker.line color. It accepts either a
                specific color or an array of numbers that are
                mapped to the colorscale relative to the max
                and min values of the array or relative to
                `marker.line.cmin` and `marker.line.cmax` if
                set.
            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
            width
                Sets the width (in px) of the lines bounding
                the marker points.
            widthsrc
                Sets the source reference on Chart Studio Cloud
                for `width`.
""",
            ),
            **kwargs,
        )
