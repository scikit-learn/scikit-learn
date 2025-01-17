import _plotly_utils.basevalidators


class Label0Validator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="label0", parent_name="pie", **kwargs):
        super(Label0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
