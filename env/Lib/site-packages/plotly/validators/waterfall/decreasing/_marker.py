import _plotly_utils.basevalidators


class MarkerValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="marker", parent_name="waterfall.decreasing", **kwargs
    ):
        super(MarkerValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Marker"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the marker color of all decreasing values.
            line
                :class:`plotly.graph_objects.waterfall.decreasi
                ng.marker.Line` instance or dict with
                compatible properties
""",
            ),
            **kwargs,
        )
