import _plotly_utils.basevalidators


class MinzoomValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="minzoom", parent_name="layout.map.layer", **kwargs):
        super(MinzoomValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            max=kwargs.pop("max", 24),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
