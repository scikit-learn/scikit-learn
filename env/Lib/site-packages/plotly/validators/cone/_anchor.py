import _plotly_utils.basevalidators


class AnchorValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="anchor", parent_name="cone", **kwargs):
        super(AnchorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["tip", "tail", "cm", "center"]),
            **kwargs,
        )
