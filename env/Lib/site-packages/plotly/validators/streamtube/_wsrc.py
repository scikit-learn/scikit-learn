import _plotly_utils.basevalidators


class WsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="wsrc", parent_name="streamtube", **kwargs):
        super(WsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
