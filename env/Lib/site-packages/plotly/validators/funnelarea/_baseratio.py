import _plotly_utils.basevalidators


class BaseratioValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="baseratio", parent_name="funnelarea", **kwargs):
        super(BaseratioValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
