import _plotly_utils.basevalidators


class AutoexpandValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="autoexpand", parent_name="layout.margin", **kwargs):
        super(AutoexpandValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
