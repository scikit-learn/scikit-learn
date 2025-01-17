import _plotly_utils.basevalidators


class ExecuteValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="execute", parent_name="layout.updatemenu.button", **kwargs
    ):
        super(ExecuteValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
