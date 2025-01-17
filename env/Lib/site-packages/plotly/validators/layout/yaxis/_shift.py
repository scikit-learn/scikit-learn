import _plotly_utils.basevalidators


class ShiftValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="shift", parent_name="layout.yaxis", **kwargs):
        super(ShiftValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
