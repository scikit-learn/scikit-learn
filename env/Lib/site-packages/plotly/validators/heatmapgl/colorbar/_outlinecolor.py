import _plotly_utils.basevalidators


class OutlinecolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self, plotly_name="outlinecolor", parent_name="heatmapgl.colorbar", **kwargs
    ):
        super(OutlinecolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
