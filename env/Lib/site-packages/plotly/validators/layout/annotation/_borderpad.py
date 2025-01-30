import _plotly_utils.basevalidators


class BorderpadValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="borderpad", parent_name="layout.annotation", **kwargs
    ):
        super(BorderpadValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+arraydraw"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
