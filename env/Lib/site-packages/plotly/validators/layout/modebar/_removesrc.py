import _plotly_utils.basevalidators


class RemovesrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="removesrc", parent_name="layout.modebar", **kwargs):
        super(RemovesrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
