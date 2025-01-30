import _plotly_utils.basevalidators


class SmoothingValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="smoothing", parent_name="contourcarpet.line", **kwargs
    ):
        super(SmoothingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            max=kwargs.pop("max", 1.3),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
