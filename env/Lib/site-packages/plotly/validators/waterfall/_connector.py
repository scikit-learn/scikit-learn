import _plotly_utils.basevalidators


class ConnectorValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="connector", parent_name="waterfall", **kwargs):
        super(ConnectorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Connector"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.waterfall.connecto
                r.Line` instance or dict with compatible
                properties
            mode
                Sets the shape of connector lines.
            visible
                Determines if connector lines are drawn.
""",
            ),
            **kwargs,
        )
