import _plotly_utils.basevalidators


class LayerdefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="layerdefaults", parent_name="layout.map", **kwargs):
        super(LayerdefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Layer"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
