import _plotly_utils.basevalidators


class TicksonValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="tickson", parent_name="layout.xaxis", **kwargs):
        super(TicksonValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            values=kwargs.pop("values", ["labels", "boundaries"]),
            **kwargs,
        )
