import _plotly_utils.basevalidators


class ShowactiveValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="showactive", parent_name="layout.updatemenu", **kwargs
    ):
        super(ShowactiveValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
