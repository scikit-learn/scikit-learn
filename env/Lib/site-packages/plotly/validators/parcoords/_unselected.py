import _plotly_utils.basevalidators


class UnselectedValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="unselected", parent_name="parcoords", **kwargs):
        super(UnselectedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Unselected"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.parcoords.unselect
                ed.Line` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
