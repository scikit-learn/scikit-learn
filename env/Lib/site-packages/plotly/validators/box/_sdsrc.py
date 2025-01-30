import _plotly_utils.basevalidators


class SdsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="sdsrc", parent_name="box", **kwargs):
        super(SdsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
