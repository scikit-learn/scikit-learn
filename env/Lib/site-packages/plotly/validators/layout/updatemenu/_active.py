import _plotly_utils.basevalidators


class ActiveValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="active", parent_name="layout.updatemenu", **kwargs):
        super(ActiveValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            min=kwargs.pop("min", -1),
            **kwargs,
        )
