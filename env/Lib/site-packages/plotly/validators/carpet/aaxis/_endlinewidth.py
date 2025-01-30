import _plotly_utils.basevalidators


class EndlinewidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="endlinewidth", parent_name="carpet.aaxis", **kwargs
    ):
        super(EndlinewidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
