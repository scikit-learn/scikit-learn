import _plotly_utils.basevalidators


class AutorangeValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="autorange", parent_name="layout.xaxis.rangeslider", **kwargs
    ):
        super(AutorangeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
