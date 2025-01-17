import _plotly_utils.basevalidators


class ScalegroupValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="scalegroup", parent_name="pie", **kwargs):
        super(ScalegroupValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
