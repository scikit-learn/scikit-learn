import _plotly_utils.basevalidators


class DepthfadeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="depthfade", parent_name="treemap.marker", **kwargs):
        super(DepthfadeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            values=kwargs.pop("values", [True, False, "reversed"]),
            **kwargs,
        )
