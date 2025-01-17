import _plotly_utils.basevalidators


class IncreasingValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="increasing", parent_name="ohlc", **kwargs):
        super(IncreasingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Increasing"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.ohlc.increasing.Li
                ne` instance or dict with compatible properties
""",
            ),
            **kwargs,
        )
