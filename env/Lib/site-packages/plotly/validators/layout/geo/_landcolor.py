import _plotly_utils.basevalidators


class LandcolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(self, plotly_name="landcolor", parent_name="layout.geo", **kwargs):
        super(LandcolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
