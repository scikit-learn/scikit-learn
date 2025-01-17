import _plotly_utils.basevalidators


class SelectedValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="selected", parent_name="scattergl", **kwargs):
        super(SelectedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Selected"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            marker
                :class:`plotly.graph_objects.scattergl.selected
                .Marker` instance or dict with compatible
                properties
            textfont
                :class:`plotly.graph_objects.scattergl.selected
                .Textfont` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
