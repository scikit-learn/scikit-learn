import _plotly_utils.basevalidators


class ParentssrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="parentssrc", parent_name="sunburst", **kwargs):
        super(ParentssrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
