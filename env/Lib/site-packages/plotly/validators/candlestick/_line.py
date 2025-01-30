import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="line", parent_name="candlestick", **kwargs):
        super(LineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Line"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            width
                Sets the width (in px) of line bounding the
                box(es). Note that this style setting can also
                be set per direction via
                `increasing.line.width` and
                `decreasing.line.width`.
""",
            ),
            **kwargs,
        )
