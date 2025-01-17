import _plotly_utils.basevalidators


class ShowticklabelsValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="showticklabels", parent_name="carpet.aaxis", **kwargs
    ):
        super(ShowticklabelsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["start", "end", "both", "none"]),
            **kwargs,
        )
