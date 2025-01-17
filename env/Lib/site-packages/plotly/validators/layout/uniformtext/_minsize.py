import _plotly_utils.basevalidators


class MinsizeValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="minsize", parent_name="layout.uniformtext", **kwargs
    ):
        super(MinsizeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
