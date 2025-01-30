import _plotly_utils.basevalidators


class YboundssrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="yboundssrc", parent_name="pointcloud", **kwargs):
        super(YboundssrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
