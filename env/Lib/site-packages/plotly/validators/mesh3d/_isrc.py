import _plotly_utils.basevalidators


class IsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="isrc", parent_name="mesh3d", **kwargs):
        super(IsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
