import _plotly_utils.basevalidators


class ColumnorderValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="columnorder", parent_name="table", **kwargs):
        super(ColumnorderValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
