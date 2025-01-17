import _plotly_utils.basevalidators


class DecreasingValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="decreasing", parent_name="candlestick", **kwargs):
        super(DecreasingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Decreasing"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fillcolor
                Sets the fill color. Defaults to a half-
                transparent variant of the line color, marker
                color, or marker line color, whichever is
                available.
            line
                :class:`plotly.graph_objects.candlestick.decrea
                sing.Line` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
