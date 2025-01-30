import _plotly_utils.basevalidators


class ArrowcolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self, plotly_name="arrowcolor", parent_name="layout.scene.annotation", **kwargs
    ):
        super(ArrowcolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
