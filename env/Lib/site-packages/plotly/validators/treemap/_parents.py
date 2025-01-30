import _plotly_utils.basevalidators


class ParentsValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="parents", parent_name="treemap", **kwargs):
        super(ParentsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
