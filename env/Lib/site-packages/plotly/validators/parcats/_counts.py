import _plotly_utils.basevalidators


class CountsValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="counts", parent_name="parcats", **kwargs):
        super(CountsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
