import _plotly_utils.basevalidators


class ViolingroupgapValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="violingroupgap", parent_name="layout", **kwargs):
        super(ViolingroupgapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
