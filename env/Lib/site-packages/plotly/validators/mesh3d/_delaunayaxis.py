import _plotly_utils.basevalidators


class DelaunayaxisValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="delaunayaxis", parent_name="mesh3d", **kwargs):
        super(DelaunayaxisValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["x", "y", "z"]),
            **kwargs,
        )
