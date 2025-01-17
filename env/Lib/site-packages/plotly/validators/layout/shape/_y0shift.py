import _plotly_utils.basevalidators


class Y0ShiftValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="y0shift", parent_name="layout.shape", **kwargs):
        super(Y0ShiftValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", -1),
            **kwargs,
        )
