import _plotly_utils.basevalidators


class StartlinewidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="startlinewidth", parent_name="carpet.aaxis", **kwargs
    ):
        super(StartlinewidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
