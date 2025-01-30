import _plotly_utils.basevalidators


class SortValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="sort", parent_name="treemap", **kwargs):
        super(SortValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
