import _plotly_utils.basevalidators


class HoleValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="hole", parent_name="layout.polar", **kwargs):
        super(HoleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
