import _plotly_utils.basevalidators


class ThresholdValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="threshold", parent_name="indicator.gauge", **kwargs
    ):
        super(ThresholdValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Threshold"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.indicator.gauge.th
                reshold.Line` instance or dict with compatible
                properties
            thickness
                Sets the thickness of the threshold line as a
                fraction of the thickness of the gauge.
            value
                Sets a treshold value drawn as a line.
""",
            ),
            **kwargs,
        )
