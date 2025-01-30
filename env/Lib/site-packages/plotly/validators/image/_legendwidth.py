import _plotly_utils.basevalidators


class LegendwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="legendwidth", parent_name="image", **kwargs):
        super(LegendwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
