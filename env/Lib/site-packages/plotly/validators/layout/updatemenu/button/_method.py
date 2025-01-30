import _plotly_utils.basevalidators


class MethodValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="method", parent_name="layout.updatemenu.button", **kwargs
    ):
        super(MethodValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            values=kwargs.pop(
                "values", ["restyle", "relayout", "animate", "update", "skip"]
            ),
            **kwargs,
        )
