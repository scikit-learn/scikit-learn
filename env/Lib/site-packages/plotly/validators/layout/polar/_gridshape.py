import _plotly_utils.basevalidators


class GridshapeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="gridshape", parent_name="layout.polar", **kwargs):
        super(GridshapeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["circular", "linear"]),
            **kwargs,
        )
