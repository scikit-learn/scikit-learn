import _plotly_utils.basevalidators


class PositionValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="position", parent_name="layout.xaxis", **kwargs):
        super(PositionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
