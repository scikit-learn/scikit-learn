import _plotly_utils.basevalidators


class XanchorValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="xanchor", parent_name="layout.shape", **kwargs):
        super(XanchorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+arraydraw"),
            **kwargs,
        )
