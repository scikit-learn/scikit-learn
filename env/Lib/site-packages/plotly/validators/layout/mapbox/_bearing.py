import _plotly_utils.basevalidators


class BearingValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="bearing", parent_name="layout.mapbox", **kwargs):
        super(BearingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
