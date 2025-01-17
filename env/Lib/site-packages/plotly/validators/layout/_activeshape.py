import _plotly_utils.basevalidators


class ActiveshapeValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="activeshape", parent_name="layout", **kwargs):
        super(ActiveshapeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Activeshape"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fillcolor
                Sets the color filling the active shape'
                interior.
            opacity
                Sets the opacity of the active shape.
""",
            ),
            **kwargs,
        )
