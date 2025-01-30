import _plotly_utils.basevalidators


class BaseValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="base", parent_name="waterfall", **kwargs):
        super(BaseValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", False),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
