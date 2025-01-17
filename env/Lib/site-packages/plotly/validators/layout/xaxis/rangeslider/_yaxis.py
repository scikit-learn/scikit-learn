import _plotly_utils.basevalidators


class YaxisValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="yaxis", parent_name="layout.xaxis.rangeslider", **kwargs
    ):
        super(YaxisValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "YAxis"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            range
                Sets the range of this axis for the
                rangeslider.
            rangemode
                Determines whether or not the range of this
                axis in the rangeslider use the same value than
                in the main plot when zooming in/out. If
                "auto", the autorange will be used. If "fixed",
                the `range` is used. If "match", the current
                range of the corresponding y-axis on the main
                subplot is used.
""",
            ),
            **kwargs,
        )
