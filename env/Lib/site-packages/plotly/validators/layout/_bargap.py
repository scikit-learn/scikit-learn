import _plotly_utils.basevalidators


class BargapValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="bargap", parent_name="layout", **kwargs):
        super(BargapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
