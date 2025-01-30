import _plotly_utils.basevalidators


class UnselectedValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="unselected", parent_name="choroplethmap", **kwargs):
        super(UnselectedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Unselected"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            marker
                :class:`plotly.graph_objects.choroplethmap.unse
                lected.Marker` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
