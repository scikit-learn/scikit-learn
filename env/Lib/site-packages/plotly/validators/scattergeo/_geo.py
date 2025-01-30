import _plotly_utils.basevalidators


class GeoValidator(_plotly_utils.basevalidators.SubplotidValidator):
    def __init__(self, plotly_name="geo", parent_name="scattergeo", **kwargs):
        super(GeoValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            dflt=kwargs.pop("dflt", "geo"),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
