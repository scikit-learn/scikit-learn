import _plotly_utils.basevalidators


class RadiusValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="radius", parent_name="layout.mapbox.layer.circle", **kwargs
    ):
        super(RadiusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
