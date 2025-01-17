import _plotly_utils.basevalidators


class GriddashValidator(_plotly_utils.basevalidators.DashValidator):
    def __init__(self, plotly_name="griddash", parent_name="layout.yaxis", **kwargs):
        super(GriddashValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            values=kwargs.pop(
                "values", ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
            ),
            **kwargs,
        )
