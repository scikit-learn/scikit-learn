import _plotly_utils.basevalidators


class TicklabelindexsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="ticklabelindexsrc", parent_name="layout.xaxis", **kwargs
    ):
        super(TicklabelindexsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
