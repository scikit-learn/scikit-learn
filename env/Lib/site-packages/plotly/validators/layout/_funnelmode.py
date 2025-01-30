import _plotly_utils.basevalidators


class FunnelmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="funnelmode", parent_name="layout", **kwargs):
        super(FunnelmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["stack", "group", "overlay"]),
            **kwargs,
        )
