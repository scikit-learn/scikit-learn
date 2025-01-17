import _plotly_utils.basevalidators


class TickwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="tickwidth", parent_name="ohlc", **kwargs):
        super(TickwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 0.5),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
