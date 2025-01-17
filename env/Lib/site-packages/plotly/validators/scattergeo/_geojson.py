import _plotly_utils.basevalidators


class GeojsonValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="geojson", parent_name="scattergeo", **kwargs):
        super(GeojsonValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
