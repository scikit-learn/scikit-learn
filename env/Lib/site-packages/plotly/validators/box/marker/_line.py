import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="line", parent_name="box.marker", **kwargs):
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
            outliercolor
                Sets the border line color of the outlier
                sample points. Defaults to marker.color
            outlierwidth
                Sets the border line width (in px) of the
                outlier sample points.
            width
                Sets the width (in px) of the lines bounding
                the marker points.
""",
            ),
            **kwargs,
        )
