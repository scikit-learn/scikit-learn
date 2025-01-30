import _plotly_utils.basevalidators


class ColumnsValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="columns", parent_name="layout.grid", **kwargs):
        super(ColumnsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
