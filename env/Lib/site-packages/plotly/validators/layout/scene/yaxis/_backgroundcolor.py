import _plotly_utils.basevalidators


class BackgroundcolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self, plotly_name="backgroundcolor", parent_name="layout.scene.yaxis", **kwargs
    ):
        super(BackgroundcolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
