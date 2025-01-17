import _plotly_utils.basevalidators


class SurfacecolorsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="surfacecolorsrc", parent_name="surface", **kwargs):
        super(SurfacecolorsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
