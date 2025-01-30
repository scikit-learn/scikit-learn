import _plotly_utils.basevalidators


class RadiusValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="radius", parent_name="densitymapbox", **kwargs):
        super(RadiusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
