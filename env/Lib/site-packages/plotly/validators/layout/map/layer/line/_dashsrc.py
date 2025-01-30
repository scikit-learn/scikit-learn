import _plotly_utils.basevalidators


class DashsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="dashsrc", parent_name="layout.map.layer.line", **kwargs
    ):
        super(DashsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
