import _plotly_utils.basevalidators


class ContoursValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="contours", parent_name="surface", **kwargs):
        super(ContoursValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Contours"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            x
                :class:`plotly.graph_objects.surface.contours.X
                ` instance or dict with compatible properties
            y
                :class:`plotly.graph_objects.surface.contours.Y
                ` instance or dict with compatible properties
            z
                :class:`plotly.graph_objects.surface.contours.Z
                ` instance or dict with compatible properties
""",
            ),
            **kwargs,
        )
