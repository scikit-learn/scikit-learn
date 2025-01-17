import _plotly_utils.basevalidators


class XsideValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="xside", parent_name="layout.grid", **kwargs):
        super(XsideValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["bottom", "bottom plot", "top plot", "top"]),
            **kwargs,
        )
