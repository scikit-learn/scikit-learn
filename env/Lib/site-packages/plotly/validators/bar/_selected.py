import _plotly_utils.basevalidators


class SelectedValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="selected", parent_name="bar", **kwargs):
        super(SelectedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Selected"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            marker
                :class:`plotly.graph_objects.bar.selected.Marke
                r` instance or dict with compatible properties
            textfont
                :class:`plotly.graph_objects.bar.selected.Textf
                ont` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
