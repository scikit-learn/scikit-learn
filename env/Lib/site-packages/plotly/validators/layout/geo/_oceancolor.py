import _plotly_utils.basevalidators


class OceancolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(self, plotly_name="oceancolor", parent_name="layout.geo", **kwargs):
        super(OceancolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
