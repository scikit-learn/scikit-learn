import _plotly_utils.basevalidators


class SourceValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="source", parent_name="sankey.link", **kwargs):
        super(SourceValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
