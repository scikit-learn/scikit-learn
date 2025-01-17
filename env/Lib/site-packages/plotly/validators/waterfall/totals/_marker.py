import _plotly_utils.basevalidators


class MarkerValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="marker", parent_name="waterfall.totals", **kwargs):
        super(MarkerValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Marker"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the marker color of all intermediate sums
                and total values.
            line
                :class:`plotly.graph_objects.waterfall.totals.m
                arker.Line` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
