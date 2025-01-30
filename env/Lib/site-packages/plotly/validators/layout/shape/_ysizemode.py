import _plotly_utils.basevalidators


class YsizemodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="ysizemode", parent_name="layout.shape", **kwargs):
        super(YsizemodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+arraydraw"),
            values=kwargs.pop("values", ["scaled", "pixel"]),
            **kwargs,
        )
