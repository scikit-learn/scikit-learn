import _plotly_utils.basevalidators


class FitboundsValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="fitbounds", parent_name="layout.geo", **kwargs):
        super(FitboundsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", [False, "locations", "geojson"]),
            **kwargs,
        )
