import _plotly_utils.basevalidators


class ActiveselectionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="activeselection", parent_name="layout", **kwargs):
        super(ActiveselectionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Activeselection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fillcolor
                Sets the color filling the active selection'
                interior.
            opacity
                Sets the opacity of the active selection.
""",
            ),
            **kwargs,
        )
