import _plotly_utils.basevalidators


class StartstandoffValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self,
        plotly_name="startstandoff",
        parent_name="layout.scene.annotation",
        **kwargs,
    ):
        super(StartstandoffValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
