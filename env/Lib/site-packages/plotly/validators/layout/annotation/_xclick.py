import _plotly_utils.basevalidators


class XclickValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="xclick", parent_name="layout.annotation", **kwargs):
        super(XclickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
