import _plotly_utils.basevalidators


class TicklabelmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="ticklabelmode", parent_name="layout.xaxis", **kwargs
    ):
        super(TicklabelmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            values=kwargs.pop("values", ["instant", "period"]),
            **kwargs,
        )
