import _plotly_utils.basevalidators


class SmithValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="smith", parent_name="layout", **kwargs):
        super(SmithValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Smith"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bgcolor
                Set the background color of the subplot
            domain
                :class:`plotly.graph_objects.layout.smith.Domai
                n` instance or dict with compatible properties
            imaginaryaxis
                :class:`plotly.graph_objects.layout.smith.Imagi
                naryaxis` instance or dict with compatible
                properties
            realaxis
                :class:`plotly.graph_objects.layout.smith.Reala
                xis` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
