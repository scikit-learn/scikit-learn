import _plotly_utils.basevalidators


class YsideValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="yside", parent_name="layout.grid", **kwargs):
        super(YsideValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["left", "left plot", "right plot", "right"]),
            **kwargs,
        )
