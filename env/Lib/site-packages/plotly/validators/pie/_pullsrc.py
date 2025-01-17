import _plotly_utils.basevalidators


class PullsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="pullsrc", parent_name="pie", **kwargs):
        super(PullsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
