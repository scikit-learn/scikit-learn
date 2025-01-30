import _plotly_utils.basevalidators


class EndlineValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="endline", parent_name="carpet.baxis", **kwargs):
        super(EndlineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
