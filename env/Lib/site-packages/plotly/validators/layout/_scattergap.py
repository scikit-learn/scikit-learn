import _plotly_utils.basevalidators


class ScattergapValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="scattergap", parent_name="layout", **kwargs):
        super(ScattergapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
