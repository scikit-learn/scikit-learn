import _plotly_utils.basevalidators


class RadiussrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="radiussrc", parent_name="densitymapbox", **kwargs):
        super(RadiussrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
