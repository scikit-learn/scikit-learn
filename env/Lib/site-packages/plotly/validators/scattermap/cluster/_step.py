import _plotly_utils.basevalidators


class StepValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="step", parent_name="scattermap.cluster", **kwargs):
        super(StepValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", -1),
            **kwargs,
        )
