import _plotly_utils.basevalidators


class YclickValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="yclick", parent_name="layout.annotation", **kwargs):
        super(YclickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
