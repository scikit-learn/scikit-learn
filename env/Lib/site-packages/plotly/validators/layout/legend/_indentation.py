import _plotly_utils.basevalidators


class IndentationValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="indentation", parent_name="layout.legend", **kwargs
    ):
        super(IndentationValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            min=kwargs.pop("min", -15),
            **kwargs,
        )
