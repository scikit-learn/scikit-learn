import _plotly_utils.basevalidators


class StartlineValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="startline", parent_name="carpet.aaxis", **kwargs):
        super(StartlineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
