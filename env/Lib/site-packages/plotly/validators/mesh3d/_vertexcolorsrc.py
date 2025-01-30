import _plotly_utils.basevalidators


class VertexcolorsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="vertexcolorsrc", parent_name="mesh3d", **kwargs):
        super(VertexcolorsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
