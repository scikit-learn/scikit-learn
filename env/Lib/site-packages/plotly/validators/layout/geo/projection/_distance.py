import _plotly_utils.basevalidators


class DistanceValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="distance", parent_name="layout.geo.projection", **kwargs
    ):
        super(DistanceValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 1.001),
            **kwargs,
        )
