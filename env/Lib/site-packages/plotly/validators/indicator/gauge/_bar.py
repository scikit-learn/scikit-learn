import _plotly_utils.basevalidators


class BarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="bar", parent_name="indicator.gauge", **kwargs):
        super(BarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Bar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the background color of the arc.
            line
                :class:`plotly.graph_objects.indicator.gauge.ba
                r.Line` instance or dict with compatible
                properties
            thickness
                Sets the thickness of the bar as a fraction of
                the total thickness of the gauge.
""",
            ),
            **kwargs,
        )
