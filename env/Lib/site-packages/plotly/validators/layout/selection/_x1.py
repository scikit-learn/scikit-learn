import _plotly_utils.basevalidators


class X1Validator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="x1", parent_name="layout.selection", **kwargs):
        super(X1Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
