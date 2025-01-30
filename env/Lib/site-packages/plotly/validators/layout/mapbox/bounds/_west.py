import _plotly_utils.basevalidators


class WestValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="west", parent_name="layout.mapbox.bounds", **kwargs
    ):
        super(WestValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
