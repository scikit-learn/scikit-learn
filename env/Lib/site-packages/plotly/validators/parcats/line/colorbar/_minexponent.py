import _plotly_utils.basevalidators


class MinexponentValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="minexponent", parent_name="parcats.line.colorbar", **kwargs
    ):
        super(MinexponentValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "colorbars"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
