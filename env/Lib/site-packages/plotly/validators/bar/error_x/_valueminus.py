import _plotly_utils.basevalidators


class ValueminusValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="valueminus", parent_name="bar.error_x", **kwargs):
        super(ValueminusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
