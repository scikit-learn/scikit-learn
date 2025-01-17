import _plotly_utils.basevalidators


class TickwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="tickwidth", parent_name="layout.yaxis", **kwargs):
        super(TickwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
