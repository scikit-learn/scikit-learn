import _plotly_utils.basevalidators


class TargetsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="targetsrc", parent_name="sankey.link", **kwargs):
        super(TargetsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
